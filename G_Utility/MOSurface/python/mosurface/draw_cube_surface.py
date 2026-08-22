"""g16_draw_cube_surface -- renders an isosurface directly from a saved
Gaussian-format .cube file, no .fchk / basis-set re-evaluation. Python
port of G_draw_cube_surface.m; see that file's docstring for the full
method notes (signed-ness detection, 'color_by' grid-matching
requirement).
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from ._common import pctile_local, read_cube_file
from ._isosurf import decimate_mesh, marching_cubes_ang
from ._render import diverging_cmap, draw_molecule_or_fallback, face_average_values, finish_3d_axes, style_kwargs

_A0 = 0.529177210544  # Bohr -> Angstrom (CODATA)

_SYM = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]


def _atomic_symbol(z):
    z = int(z)
    return _SYM[z - 1] if 1 <= z <= len(_SYM) else f"Z{z}"


def g16_draw_cube_surface(cubefile, color_by="", isovalue=None, surface_style="transparent",
                           face_alpha=None, pos_color=(0.10, 0.40, 0.85), neg_color=(0.85, 0.15, 0.10),
                           max_vertices=5000, show_molecule=True, show_labels=False, atom_scale=0.35,
                           title="", ax=None, verbose=True):
    """Renders an isosurface directly from a .cube file (as saved by any
    of the three g16_draw_*_surface functions' save_cube option, or by
    Gaussian's own cubegen).

    color_by : path to a SECOND .cube file. If given, the isosurface
        SHAPE still comes from cubefile, but coloured by trilinearly
        interpolating the second cube's field -- GaussView's classic
        "ESP mapped on electron density" from two independently-saved
        cube files. The two cubes must share an identical grid.
    isovalue : isosurface level. Default: guessed from cubefile's title
        (0.001 density-like, 0.02 MO-cube); no default for anything else.

    Returns a dict h with keys 'pos'/'neg' (Poly3DCollection or None).
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    face_alpha_explicit = face_alpha is not None
    if face_alpha is None:
        face_alpha = 1.0 if surface_style == "solid" else 0.55
    if face_alpha_explicit and surface_style == "grid" and verbose:
        print("g16_draw_cube_surface: face_alpha has no effect with surface_style='grid'.")

    color_mode = bool(color_by)

    cube = read_cube_file(cubefile)
    gx, gy, gz, V = cube["gx"], cube["gy"], cube["gz"], cube["V"]
    atomic_numbers, xyz_bohr = cube["atomic_numbers"], cube["xyz_bohr"]
    orbital_indices, title_line = cube["orbital_indices"], cube["title_line"]
    is_mo_cube = orbital_indices is not None

    if verbose:
        tag = f" (MO cube, orbital {list(orbital_indices)})" if is_mo_cube else ""
        print(f"g16_draw_cube_surface: {cubefile} -- grid {len(gx)} x {len(gy)} x {len(gz)}, "
              f"{len(atomic_numbers)} atoms{tag}")

    if isovalue is not None:
        isoval = abs(isovalue)
    elif is_mo_cube:
        isoval = 0.02
        if verbose:
            print("  isovalue not given -- guessed 0.02 (MO cube).")
    elif "density" in title_line.lower():
        isoval = 0.001
        if verbose:
            print(f"  isovalue not given -- guessed 0.001 (density-like title: {title_line!r}).")
    else:
        raise ValueError(
            f"g16_draw_cube_surface: isovalue was not given, and the field type could not be "
            f"guessed from the cube's title ({title_line!r}). Pass isovalue explicitly."
        )

    if not color_mode and ("electrostatic potential" in title_line.lower() or "esp " in title_line.lower()):
        print(f"g16_draw_cube_surface: WARNING {cubefile} appears to be an ESP field (title: {title_line!r}), "
              f"isosurfaced here directly at isovalue={isoval:.4g} -- this is a surface of CONSTANT ESP "
              "VALUE, not \"ESP mapped on the density envelope\" (use color_by for that).")

    if color_mode:
        signed_field = False
    elif is_mo_cube:
        signed_field = True
    else:
        vtol = 1e-6 * np.max(np.abs(V))
        signed_field = np.any(V > vtol) and np.any(V < -vtol)

    Xa, Ya, Za = gx * _A0, gy * _A0, gz * _A0

    verts_pos, faces_pos, normals_pos = marching_cubes_ang(Xa, Ya, Za, V, isoval)
    if signed_field:
        verts_neg, faces_neg, normals_neg = marching_cubes_ang(Xa, Ya, Za, V, -isoval)
    else:
        verts_neg, faces_neg, normals_neg = None, None, None

    if verts_pos is None and (not signed_field or verts_neg is None):
        raise ValueError(
            f"g16_draw_cube_surface: no isosurface found at isovalue={isoval:.4g} -- the cube's value "
            f"range is [{V.min():.4g}, {V.max():.4g}]."
        )

    n0 = verts_pos.shape[0] if verts_pos is not None else 0
    if verts_pos is not None and n0 > max_vertices:
        dec = decimate_mesh(verts_pos, faces_pos, max_vertices, normals=normals_pos)
        verts_pos, faces_pos, normals_pos = dec["vertices"], dec["faces"], dec["normals"]
        if verbose:
            print(f"  Positive lobe decimated: {n0} -> {verts_pos.shape[0]} vertices")
    if signed_field and verts_neg is not None and verts_neg.shape[0] > max_vertices:
        n0n = verts_neg.shape[0]
        dec = decimate_mesh(verts_neg, faces_neg, max_vertices, normals=normals_neg)
        verts_neg, faces_neg, normals_neg = dec["vertices"], dec["faces"], dec["normals"]
        if verbose:
            print(f"  Negative lobe decimated: {n0n} -> {verts_neg.shape[0]} vertices")

    cdata = None
    if color_mode:
        cube2 = read_cube_file(color_by)
        gx2, gy2, gz2, V2 = cube2["gx"], cube2["gy"], cube2["gz"], cube2["V"]
        same_grid = (len(gx2) == len(gx) and len(gy2) == len(gy) and len(gz2) == len(gz) and
                     np.max(np.abs(gx2 - gx)) < 1e-4 and np.max(np.abs(gy2 - gy)) < 1e-4 and
                     np.max(np.abs(gz2 - gz)) < 1e-4)
        if not same_grid:
            raise ValueError(
                f"g16_draw_cube_surface: color_by ({color_by}) does not share an identical grid with "
                f"{cubefile} (origin/spacing/point-counts must match exactly)."
            )
        # V2 has meshgrid convention V2[iy,ix,iz]; RegularGridInterpolator wants
        # values indexed (x,y,z) matching its own axis order.
        # bounds_error=False, fill_value=None extrapolates (linearly) rather
        # than returning NaN for the rare marching-cubes vertex that lands
        # exactly on -- or a floating-point hair outside -- the last grid
        # cell edge; since such a vertex is by construction at/adjacent to
        # the domain boundary, linear extrapolation there is effectively
        # identical to the true boundary value (unlike MATLAB's interp3,
        # which would leave a handful of such vertices as NaN too).
        interp = RegularGridInterpolator((gx, gy, gz), np.transpose(V2, (1, 0, 2)),
                                          bounds_error=False, fill_value=None)
        verts_bohr = verts_pos / _A0
        cdata = interp(verts_bohr)
        if verbose:
            print(f"  color_by field interpolated at {cdata.size} vertices: "
                  f"range [{np.nanmin(cdata):.4g}, {np.nanmax(cdata):.4g}]"
                  + (f" ({np.isnan(cdata).sum()} NaN)" if np.any(np.isnan(cdata)) else ""))

    if ax is None:
        fig = plt.figure()
        fig.canvas.manager.set_window_title("Cube surface")
        ax = fig.add_subplot(111, projection="3d")

    h = {"pos": None, "neg": None}

    if color_mode:
        cmap = diverging_cmap(neg_color, pos_color)
        vmax = pctile_local(np.abs(cdata), 99)
        if vmax == 0:
            vmax = 1.0
        face_vals = face_average_values(cdata, faces_pos)
        face_colors = cmap((face_vals + vmax) / (2 * vmax))
        if surface_style == "grid":
            h["pos"] = Poly3DCollection(verts_pos[faces_pos], facecolors=(0, 0, 0, 0),
                                         edgecolors=face_colors, linewidths=0.6, alpha=1.0)
        else:
            h["pos"] = Poly3DCollection(verts_pos[faces_pos], facecolors=face_colors,
                                         edgecolors=(0, 0, 0, 0), alpha=face_alpha, shade=False)
        ax.add_collection3d(h["pos"])

        import matplotlib.cm as mcm
        import matplotlib.colors as mcolors
        sm = mcm.ScalarMappable(norm=mcolors.Normalize(vmin=-vmax, vmax=vmax), cmap=cmap)
        cb = plt.colorbar(sm, ax=ax)
        cb.set_label("color_by field value")
    else:
        if verts_pos is not None:
            h["pos"] = Poly3DCollection(verts_pos[faces_pos], **style_kwargs(surface_style, pos_color, face_alpha))
            ax.add_collection3d(h["pos"])
        if signed_field and verts_neg is not None:
            h["neg"] = Poly3DCollection(verts_neg[faces_neg], **style_kwargs(surface_style, neg_color, face_alpha))
            ax.add_collection3d(h["neg"])

    finish_3d_axes(ax, xyz_bohr * _A0, verts_pos, verts_neg)

    if show_molecule:
        mol = type("Mol", (), {})()
        mol.Natoms = len(atomic_numbers)
        mol.symbols = [_atomic_symbol(z) for z in atomic_numbers]
        mol.xyz = xyz_bohr * _A0
        mol.filename = cubefile
        draw_molecule_or_fallback(ax, mol, show_labels, atom_scale, "g16_draw_cube_surface")

    if not title:
        title = f"{title_line} (iso={isoval:.3g})"
    ax.set_title(title, fontsize=11)

    return h
