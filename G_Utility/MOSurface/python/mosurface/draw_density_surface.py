"""g16_draw_density_surface -- total/difference/spin electron density as a
real-space isosurface or 2D contour. Python port of
G_draw_density_surface.m; see that file's docstring for the full method
notes.
"""

import time

import numpy as np

from ._common import (check_supported_shells, is_clipped, pctile_local,
                       read_aobasis_from_fchk, read_fchk_scalar_or_section,
                       write_cube_file)
from ._eval import eval_density_on_grid
from ._isosurf import marching_cubes_ang
from ._plot2d import draw_2d_atoms, plane_basis
from ._render import diverging_cmap, draw_molecule_or_fallback, finish_3d_axes, sequential_cmap, style_kwargs

_A0 = 0.529177210544  # Bohr -> Angstrom (CODATA)


def _validate_and_prepare(data, argname):
    for req in ("filename", "mol", "Nbasis", "Nbasis_indep", "alpha_MO_coeff", "xyz_bohr", "Nalpha", "Nbeta"):
        if not hasattr(data, req):
            raise ValueError(f"g16_draw_density_surface: {argname} must be the Struct returned by g16_fchk_read (missing '{req}').")
    if data.Nalpha != data.Nbeta:
        raise ValueError(
            f"g16_draw_density_surface: {argname} is an open-shell calculation "
            f"(Nalpha={data.Nalpha}, Nbeta={data.Nbeta}) -- only closed-shell is supported here."
        )
    aobasis = read_aobasis_from_fchk(data.filename)
    if aobasis["shell_types"].size == 0:
        raise ValueError(f"g16_draw_density_surface: no basis-set sections found in {data.filename} ({argname}).")
    check_supported_shells(aobasis["shell_types"])
    alpha_MO = np.asarray(data.alpha_MO_coeff).reshape(data.Nbasis, data.Nbasis_indep, order="F")
    occ_coeff = alpha_MO[:, :data.Nalpha]
    return aobasis, occ_coeff


def _validate_and_prepare_spin(data, argname):
    for req in ("filename", "mol", "Nbasis", "Nbasis_indep", "alpha_MO_coeff", "xyz_bohr", "Nalpha", "Nbeta"):
        if not hasattr(data, req):
            raise ValueError(f"g16_draw_density_surface: {argname} must be the Struct returned by g16_fchk_read (missing '{req}').")
    aobasis = read_aobasis_from_fchk(data.filename)
    if aobasis["shell_types"].size == 0:
        raise ValueError(f"g16_draw_density_surface: no basis-set sections found in {data.filename} ({argname}).")
    check_supported_shells(aobasis["shell_types"])
    alpha_MO = np.asarray(data.alpha_MO_coeff).reshape(data.Nbasis, data.Nbasis_indep, order="F")
    occ_alpha = alpha_MO[:, :data.Nalpha]

    beta_vec = read_fchk_scalar_or_section(data.filename, "Beta MO coefficients")
    if beta_vec is None:
        raise ValueError(
            f"g16_draw_density_surface: spin_density requires an unrestricted (UHF/UKS) "
            f"calculation, but no 'Beta MO coefficients' section was found in {data.filename} ({argname})."
        )
    expected = data.Nbasis * data.Nbasis_indep
    if beta_vec.size != expected:
        raise ValueError(
            f"g16_draw_density_surface: 'Beta MO coefficients' section in {data.filename} has "
            f"{beta_vec.size} entries, expected {expected} (Nbasis*Nbasis_indep)."
        )
    beta_MO = beta_vec.reshape(data.Nbasis, data.Nbasis_indep, order="F")
    occ_beta = beta_MO[:, :data.Nbeta]
    return aobasis, occ_alpha, occ_beta


def g16_draw_density_surface(data, compare_to=None, spin_density=False, mode="surface",
                              plane="auto", plane_offset=0.0, isovalue=0.001,
                              grid_spacing=0.15, padding=4.0, save_cube="",
                              pos_color=(0.10, 0.40, 0.85), neg_color=(0.85, 0.15, 0.10),
                              face_alpha=0.55, surface_style="transparent",
                              show_molecule=True, show_labels=False, atom_scale=0.35,
                              title="", ax=None, verbose=True, padding_explicit=False):
    """Renders the total electron density (or, with compare_to, a density
    DIFFERENCE; or, with spin_density, the SPIN density alpha-minus-beta)
    as a real-space isosurface or 2D contour.

    See G_draw_density_surface.m's docstring for the full parameter list
    -- names here are the same options in snake_case.

    Returns a dict h with keys 'pos'/'neg' (Poly3DCollection or None) in
    'surface' mode, 'contour'/'zero' (QuadContourSet/None) in 'contour'
    mode.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    compare_mode = compare_to is not None
    signed_mode = compare_mode or spin_density
    mode = mode.lower()
    plane_choice = plane.lower()
    isoval = abs(isovalue)
    spacing = grid_spacing
    pad = padding
    face_alpha_explicit = face_alpha != 0.55

    if surface_style == "solid" and not face_alpha_explicit:
        face_alpha = 1.0
    if face_alpha_explicit and surface_style == "grid" and verbose:
        print("g16_draw_density_surface: 'face_alpha' has no effect with surface_style='grid'.")
    if save_cube and mode == "contour":
        raise ValueError("g16_draw_density_surface: save_cube needs the full 3D grid built in mode='surface'.")
    if spin_density and compare_mode:
        raise ValueError("g16_draw_density_surface: spin_density and compare_to are mutually exclusive.")

    if spin_density:
        aobasis1, occ_alpha, occ_beta = _validate_and_prepare_spin(data, "data")
    else:
        aobasis1, occ1 = _validate_and_prepare(data, "data")
        if compare_mode:
            aobasis2, occ2 = _validate_and_prepare(compare_to, "compare_to")

    xyz_bohr = np.asarray(data.xyz_bohr)

    if ax is None:
        fig = plt.figure()
        fig.canvas.manager.set_window_title("Density surface")
        ax = fig.add_subplot(111, projection="3d" if mode == "surface" else None)

    h = {"pos": None, "neg": None, "contour": None, "zero": None}

    def eval_rho(pts):
        if spin_density:
            return (eval_density_on_grid(aobasis1, occ_alpha, pts, 1)
                    - eval_density_on_grid(aobasis1, occ_beta, pts, 1))
        rho = eval_density_on_grid(aobasis1, occ1, pts)
        if compare_mode:
            rho = eval_density_on_grid(aobasis2, occ2, pts) - rho
        return rho

    if mode == "surface":
        max_grid, max_attempts, pad_growth = 30_000_000, 4, 1.6
        attempt = 0
        clipped = False
        while True:
            attempt += 1
            lo = xyz_bohr.min(axis=0) - pad
            hi = xyz_bohr.max(axis=0) + pad
            gx = np.arange(lo[0], hi[0] + spacing / 2, spacing)
            gy = np.arange(lo[1], hi[1] + spacing / 2, spacing)
            gz = np.arange(lo[2], hi[2] + spacing / 2, spacing)
            Ngrid = len(gx) * len(gy) * len(gz)

            if Ngrid > max_grid:
                if attempt == 1:
                    raise ValueError(f"g16_draw_density_surface: requested grid has {Ngrid} points, too large.")
                print(f"g16_draw_density_surface: stopped auto-expanding padding at {pad:.1f} Bohr.")
                break

            if verbose:
                print(f"g16_draw_density_surface: grid {len(gx)} x {len(gy)} x {len(gz)} = {Ngrid} points "
                      f"(spacing={spacing:.2f} Bohr, padding={pad:.1f} Bohr)")

            X, Y, Z = np.meshgrid(gx, gy, gz)
            grid_pts_bohr = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

            t0 = time.time()
            V = eval_rho(grid_pts_bohr).reshape(X.shape)
            if verbose:
                print(f"  Density evaluation: {time.time()-t0:.1f} s")

            Xa, Ya, Za = gx * _A0, gy * _A0, gz * _A0
            verts_pos, faces_pos, _ = marching_cubes_ang(Xa, Ya, Za, V, isoval)
            if signed_mode:
                verts_neg, faces_neg, _ = marching_cubes_ang(Xa, Ya, Za, V, -isoval)
            else:
                verts_neg, faces_neg = None, None

            clipped = (is_clipped(verts_pos, Xa, Ya, Za, 1.5 * spacing * _A0) or
                       is_clipped(verts_neg, Xa, Ya, Za, 1.5 * spacing * _A0))

            if not clipped or padding_explicit or attempt >= max_attempts:
                break
            if verbose:
                print(f"  Isosurface touches the grid boundary -- expanding padding {pad:.1f} -> {pad*pad_growth:.1f} Bohr.")
            pad *= pad_growth

        if clipped:
            print(f"g16_draw_density_surface: WARNING isosurface appears clipped (padding={pad:.1f} Bohr).")

        if save_cube:
            if spin_density:
                title_line = "Spin density (alpha - beta) -- G_Utility export"
            elif compare_mode:
                title_line = "Density difference (compare_to - data) -- G_Utility export"
            else:
                title_line = "Total electron density -- G_Utility export"
            comment_line = f"rho(r), electrons/Bohr^3, grid spacing={spacing:.3f} Bohr"
            write_cube_file(save_cube, title_line, comment_line, data.mol.Z, xyz_bohr, gx, gy, gz, V)
            if verbose:
                print(f"  Cube file written: {save_cube}")

        if verts_pos is not None:
            h["pos"] = Poly3DCollection(verts_pos[faces_pos], **style_kwargs(surface_style, pos_color, face_alpha))
            ax.add_collection3d(h["pos"])
        if signed_mode and verts_neg is not None:
            h["neg"] = Poly3DCollection(verts_neg[faces_neg], **style_kwargs(surface_style, neg_color, face_alpha))
            ax.add_collection3d(h["neg"])
        if verts_pos is None and (not signed_mode or verts_neg is None):
            print(f"g16_draw_density_surface: WARNING no isosurface found at isovalue={isoval:.4g}.")

        finish_3d_axes(ax, xyz_bohr * _A0, verts_pos, verts_neg)
        if show_molecule:
            draw_molecule_or_fallback(ax, data.mol, show_labels, atom_scale, "g16_draw_density_surface")

    else:  # contour
        center, u, v, n = plane_basis(xyz_bohr, plane_choice)
        offset_bohr = plane_offset / _A0
        atom_uv = (xyz_bohr - center) @ np.column_stack([u, v])
        lo_uv = atom_uv.min(axis=0) - pad
        hi_uv = atom_uv.max(axis=0) + pad
        gs = np.arange(lo_uv[0], hi_uv[0] + spacing / 2, spacing)
        gt = np.arange(lo_uv[1], hi_uv[1] + spacing / 2, spacing)
        Ngrid = len(gs) * len(gt)
        if Ngrid > 4_000_000:
            raise ValueError(f"g16_draw_density_surface: requested contour grid has {Ngrid} points, too large.")

        if verbose:
            print(f"g16_draw_density_surface: contour grid {len(gs)} x {len(gt)} = {Ngrid} points "
                  f"(spacing={spacing:.2f} Bohr, padding={pad:.1f} Bohr, plane={plane_choice}, offset={plane_offset:.2f} A)")

        S, T = np.meshgrid(gs, gt)
        grid_pts_bohr = (center + offset_bohr * n) + S.ravel()[:, None] * u + T.ravel()[:, None] * v

        t0 = time.time()
        V2 = eval_rho(grid_pts_bohr).reshape(S.shape)
        if verbose:
            print(f"  Density evaluation: {time.time()-t0:.1f} s")

        Sa, Ta = S * _A0, T * _A0

        if signed_mode:
            cmap = diverging_cmap(neg_color, pos_color)
            vmax = pctile_local(np.abs(V2), 98)
            if vmax == 0:
                vmax = 1.0
            clow, chigh = -vmax, vmax
        else:
            cmap = sequential_cmap(pos_color)
            vmax = pctile_local(V2, 98)
            if vmax <= 0:
                vmax = 1.0
            clow, chigh = 0.0, vmax

        levels = np.linspace(clow, chigh, 21)
        h["contour"] = ax.contourf(Sa, Ta, V2, levels=levels, cmap=cmap, vmin=clow, vmax=chigh, extend="both")

        if signed_mode and vmax > 0 and np.any(V2 > 0) and np.any(V2 < 0):
            h["zero"] = ax.contour(Sa, Ta, V2, levels=[0], colors=[(0.15, 0.15, 0.15)], linewidths=1.3)

        cb = plt.colorbar(h["contour"], ax=ax)
        if spin_density:
            cb.set_label("Spin density, alpha - beta (e/Bohr^3)")
        elif compare_mode:
            cb.set_label("Density difference (e/Bohr^3)")
        else:
            cb.set_label("Electron density (e/Bohr^3)")

        ax.set_aspect("equal")
        ax.set_xlabel("in-plane u (A)")
        ax.set_ylabel("in-plane v (A)")

        if show_molecule:
            draw_2d_atoms(ax, data.mol.symbols, atom_uv * _A0, xyz_bohr * _A0)

    if not title:
        if spin_density:
            base_str = "Spin density ($\\alpha$-$\\beta$)"
        elif compare_mode:
            base_str = "Density difference"
        else:
            base_str = "Total density"
        if mode == "surface":
            title = f"{base_str}, iso={'$\\pm$' if signed_mode else ''}{isoval:.3g}"
        else:
            title = f"{base_str}, plane={plane_choice}"
    ax.set_title(title, fontsize=11)

    return h
