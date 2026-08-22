"""g16_draw_mo_surface -- real-space molecular-orbital isosurface (or 2D
nodal contour map), evaluated from the raw .fchk basis-set data, overlaid
on the CPK molecule. Python port of G_draw_mo_surface.m; see that file's
docstring for the full method/normalization/validation notes.
"""

import re
import time

import numpy as np

from ._common import check_supported_shells, is_clipped, read_aobasis_from_fchk, write_cube_file
from ._eval import eval_mo_on_grid
from ._isosurf import marching_cubes_ang
from ._plot2d import draw_2d_atoms, plane_basis
from ._render import diverging_cmap, draw_molecule_or_fallback, finish_3d_axes, style_kwargs

_A0 = 0.529177210544  # Bohr -> Angstrom (CODATA)


def _resolve_mo_index(mo_spec, homo_idx, n_basis):
    if isinstance(mo_spec, (int, np.integer)):
        idx = int(mo_spec)
    elif isinstance(mo_spec, str):
        s = mo_spec.strip().upper()
        m = re.match(r'^(HOMO|LUMO)([+-]\d+)?$', s)
        if not m:
            raise ValueError("resolve_mo_index: mo_index string must be numeric, or 'HOMO'/'LUMO'/'HOMO-n'/'LUMO+n'.")
        if homo_idx is None:
            raise ValueError(f"resolve_mo_index: mo_index={s!r} requires data.HOMO_idx, which is missing.")
        base, offs = m.group(1), m.group(2)
        n = int(offs) if offs else 0
        idx = homo_idx + n if base == "HOMO" else homo_idx + 1 + n
    else:
        raise TypeError("resolve_mo_index: mo_index must be an int or a string like 'HOMO'/'LUMO'/'HOMO-1'/'LUMO+2'.")
    if idx < 1 or idx > n_basis:
        raise ValueError(f"resolve_mo_index: resolved MO index {idx} is out of range [1, {n_basis}].")
    return idx


def g16_draw_mo_surface(data, mo_index, mode="surface", plane="auto", plane_offset=0.0,
                         isovalue=0.02, grid_spacing=0.15, padding=4.0, save_cube="",
                         pos_color=(0.10, 0.40, 0.85), neg_color=(0.85, 0.15, 0.10),
                         face_alpha=0.55, surface_style="transparent",
                         show_molecule=True, show_labels=False, atom_scale=0.35,
                         title="", ax=None, verbose=True, padding_explicit=None):
    """Renders a molecular-orbital real-space isosurface (positive/
    negative lobes) or 2D nodal contour map.

    data : Struct/dict as returned by G16parser's g16_fchk_read. Required
        attributes: filename, mol, Nbasis, alpha_MO_coeff, xyz_bohr.
        Optional: HOMO_idx, alpha_orb_energies (title/HOMO-LUMO tagging).
    mo_index : positive int (1-based MO column), or 'HOMO'/'LUMO'/
        'HOMO-n'/'LUMO+n' (resolved against data.HOMO_idx).

    See G_draw_mo_surface.m's docstring for the full parameter list --
    names here are the same options in snake_case. padding_explicit: set
    True/False explicitly if calling programmatically to control the
    auto-padding-expansion behaviour; defaults to None, meaning "explicit
    iff padding was passed a non-default value" is not tracked here (pass
    it explicitly if you need MATLAB's exact auto-detection semantics).

    Returns a dict h with keys 'pos'/'neg' (Poly3DCollection or None) in
    'surface' mode, 'contour'/'zero' (QuadContourSet/None) in 'contour'
    mode.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    mode = mode.lower()
    plane_choice = plane.lower()
    isoval = abs(isovalue)
    spacing = grid_spacing
    pad = padding
    if padding_explicit is None:
        padding_explicit = False
    face_alpha_explicit = face_alpha != 0.55

    if surface_style == "solid" and not face_alpha_explicit:
        face_alpha = 1.0
    if face_alpha_explicit and surface_style == "grid" and verbose:
        print("g16_draw_mo_surface: 'face_alpha' has no effect with surface_style='grid' (no face is drawn, only mesh edges).")

    if save_cube and mode == "contour":
        raise ValueError("g16_draw_mo_surface: save_cube needs the full 3D grid built in mode='surface' -- not available in mode='contour'.")

    for req in ("filename", "mol", "Nbasis", "alpha_MO_coeff", "xyz_bohr"):
        if not hasattr(data, req):
            raise ValueError(f"g16_draw_mo_surface: data must be the Struct returned by g16_fchk_read (missing '{req}').")

    aobasis = read_aobasis_from_fchk(data.filename)
    if aobasis["shell_types"].size == 0:
        raise ValueError(f"g16_draw_mo_surface: no basis-set sections found in {data.filename}.")
    check_supported_shells(aobasis["shell_types"])

    homo_idx = getattr(data, "HOMO_idx", None)
    idx = _resolve_mo_index(mo_index, homo_idx, data.Nbasis)

    alpha_MO = np.asarray(data.alpha_MO_coeff).reshape(data.Nbasis, data.Nbasis, order="F")
    mo_col = alpha_MO[:, idx - 1:idx]

    xyz_bohr = np.asarray(data.xyz_bohr)

    if ax is None:
        fig = plt.figure()
        fig.canvas.manager.set_window_title("MO surface")
        ax = fig.add_subplot(111, projection="3d" if mode == "surface" else None)

    h = {"pos": None, "neg": None, "contour": None, "zero": None}

    if mode == "surface":
        max_grid = 30_000_000
        max_attempts = 4
        pad_growth = 1.6

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
                    raise ValueError(
                        f"g16_draw_mo_surface: requested grid has {Ngrid} points "
                        f"({len(gx)} x {len(gy)} x {len(gz)}), too large. Increase "
                        "grid_spacing or decrease padding."
                    )
                print(f"g16_draw_mo_surface: stopped auto-expanding padding at {pad:.1f} Bohr to stay under the {max_grid}-point grid limit.")
                break

            if verbose:
                print(f"g16_draw_mo_surface: MO {idx}, grid {len(gx)} x {len(gy)} x {len(gz)} "
                      f"= {Ngrid} points (spacing={spacing:.2f} Bohr, padding={pad:.1f} Bohr)")

            X, Y, Z = np.meshgrid(gx, gy, gz)
            grid_pts_bohr = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

            t0 = time.time()
            mo_vals = eval_mo_on_grid(aobasis, mo_col, grid_pts_bohr)
            V = mo_vals.reshape(X.shape)
            if verbose:
                print(f"  MO evaluation: {time.time()-t0:.1f} s")

            Xa, Ya, Za = gx * _A0, gy * _A0, gz * _A0
            verts_pos, faces_pos, normals_pos = marching_cubes_ang(Xa, Ya, Za, V, isoval)
            verts_neg, faces_neg, normals_neg = marching_cubes_ang(Xa, Ya, Za, V, -isoval)

            clipped = (is_clipped(verts_pos, Xa, Ya, Za, 1.5 * spacing * _A0) or
                       is_clipped(verts_neg, Xa, Ya, Za, 1.5 * spacing * _A0))

            if not clipped or padding_explicit or attempt >= max_attempts:
                break
            if verbose:
                print(f"  Isosurface touches the grid boundary -- expanding padding {pad:.1f} -> {pad*pad_growth:.1f} Bohr and retrying.")
            pad *= pad_growth

        if clipped:
            if padding_explicit:
                print(f"g16_draw_mo_surface: WARNING the isosurface appears clipped by the grid boundary (padding={pad:.1f} Bohr). Try a larger padding.")
            else:
                print(f"g16_draw_mo_surface: WARNING the isosurface still appears clipped after auto-expanding padding to {pad:.1f} Bohr ({attempt} attempts).")

        if save_cube:
            title_line = f"MO {idx} amplitude -- G_Utility export"
            comment_line = f"psi(r), 1/Bohr^1.5, grid spacing={spacing:.3f} Bohr"
            write_cube_file(save_cube, title_line, comment_line, data.mol.Z, xyz_bohr, gx, gy, gz, V, [idx])
            if verbose:
                print(f"  Cube file written: {save_cube}")

        if verts_pos is not None:
            h["pos"] = Poly3DCollection(verts_pos[faces_pos], **style_kwargs(surface_style, pos_color, face_alpha))
            ax.add_collection3d(h["pos"])
        if verts_neg is not None:
            h["neg"] = Poly3DCollection(verts_neg[faces_neg], **style_kwargs(surface_style, neg_color, face_alpha))
            ax.add_collection3d(h["neg"])
        if verts_pos is None and verts_neg is None:
            print(f"g16_draw_mo_surface: WARNING no isosurface found at isovalue={isoval:.4g}.")

        finish_3d_axes(ax, xyz_bohr * _A0, verts_pos, verts_neg)

        if show_molecule:
            draw_molecule_or_fallback(ax, data.mol, show_labels, atom_scale)

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
            raise ValueError(f"g16_draw_mo_surface: requested contour grid has {Ngrid} points, too large.")

        if verbose:
            print(f"g16_draw_mo_surface: MO {idx}, contour grid {len(gs)} x {len(gt)} = {Ngrid} points "
                  f"(spacing={spacing:.2f} Bohr, padding={pad:.1f} Bohr, plane={plane_choice}, offset={plane_offset:.2f} A)")

        S, T = np.meshgrid(gs, gt)
        grid_pts_bohr = (center + offset_bohr * n) + S.ravel()[:, None] * u + T.ravel()[:, None] * v

        t0 = time.time()
        mo_vals = eval_mo_on_grid(aobasis, mo_col, grid_pts_bohr)
        V2 = mo_vals.reshape(S.shape)
        if verbose:
            print(f"  MO evaluation: {time.time()-t0:.1f} s")

        vmax = np.max(np.abs(V2))
        if vmax < 1e-3:
            print(f"g16_draw_mo_surface: WARNING the MO amplitude is essentially zero (max {vmax:.2e} a.u.) "
                  f"everywhere on this plane (offset={plane_offset:.2f} A). Try a nonzero plane_offset.")

        Sa, Ta = S * _A0, T * _A0

        cmap = diverging_cmap(neg_color, pos_color)
        if vmax == 0:
            vmax = 1.0
        levels = np.linspace(-vmax, vmax, 21)
        h["contour"] = ax.contourf(Sa, Ta, V2, levels=levels, cmap=cmap, vmin=-vmax, vmax=vmax)

        if vmax > 0 and np.any(V2 > 0) and np.any(V2 < 0):
            h["zero"] = ax.contour(Sa, Ta, V2, levels=[0], colors=[(0.15, 0.15, 0.15)], linewidths=1.3)

        cb = plt.colorbar(h["contour"], ax=ax)
        cb.set_label("MO amplitude (a.u.)")

        ax.set_aspect("equal")
        ax.set_xlabel("in-plane u (A)")
        ax.set_ylabel("in-plane v (A)")

        if show_molecule:
            draw_2d_atoms(ax, data.mol.symbols, atom_uv * _A0, xyz_bohr * _A0)

    if not title:
        tag = ""
        if homo_idx is not None:
            if idx == homo_idx:
                tag = " (HOMO)"
            elif idx == homo_idx + 1:
                tag = " (LUMO)"
            elif idx < homo_idx:
                tag = f" (HOMO-{homo_idx - idx})"
            else:
                tag = f" (LUMO+{idx - homo_idx - 1})"
        e_str = ""
        orb_e = getattr(data, "alpha_orb_energies", None)
        if orb_e is not None and idx <= len(orb_e):
            e_str = f", E = {orb_e[idx-1]:.4f} Ha"
        if mode == "surface":
            title = f"MO {idx}{tag}{e_str}, iso=$\\pm${isoval:.3g}"
        else:
            title = f"MO {idx}{tag}{e_str}, plane={plane_choice}"
    ax.set_title(title, fontsize=11)

    return h
