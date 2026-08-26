"""g16_draw_esp_surface -- electrostatic potential mapped onto the total
electron density isosurface. Python port of G_draw_esp_surface.m; see
that file's docstring for the full method/accuracy notes (analytic
McMurchie-Davidson vs. numerical Coulomb-grid-sum fallback).
"""

import time

import numpy as np

from ._common import (check_supported_shells, is_clipped, pctile_local,
                       read_aobasis_from_fchk, read_fchk_scalar_or_section,
                       write_cube_file)
from ._eval import eval_density_on_grid, eval_esp_analytic
from ._isosurf import decimate_mesh, marching_cubes_ang
from ._render import diverging_cmap, draw_molecule_or_fallback, face_average_values, finish_3d_axes

_A0 = 0.529177210544  # Bohr -> Angstrom (CODATA)


def _read_nuclear_charges(filename, natoms):
    Z = read_fchk_scalar_or_section(filename, "Nuclear charges")
    if Z is None:
        raise ValueError(f"g16_draw_esp_surface: 'Nuclear charges' section not found in {filename}.")
    if Z.size != natoms:
        raise ValueError(
            f"g16_draw_esp_surface: 'Nuclear charges' section in {filename} has {Z.size} "
            f"entries, expected {natoms} (Natoms)."
        )
    return Z


def _validate_and_prepare(data, argname):
    for req in ("filename", "mol", "Nbasis", "Nbasis_indep", "alpha_MO_coeff", "xyz_bohr", "Nalpha", "Nbeta", "Nelec"):
        if not hasattr(data, req):
            raise ValueError(f"g16_draw_esp_surface: {argname} must be the Struct returned by g16_fchk_read (missing '{req}').")
    if data.Nalpha != data.Nbeta:
        raise ValueError(f"g16_draw_esp_surface: {argname} is an open-shell calculation -- only closed-shell is supported.")

    aobasis = read_aobasis_from_fchk(data.filename)
    if aobasis["shell_types"].size == 0:
        raise ValueError(f"g16_draw_esp_surface: no basis-set sections found in {data.filename} ({argname}).")
    check_supported_shells(aobasis["shell_types"])

    alpha_MO = np.asarray(data.alpha_MO_coeff).reshape(data.Nbasis, data.Nbasis_indep, order="F")
    occ_coeff = alpha_MO[:, :data.Nalpha]
    P_density = 2 * (occ_coeff @ occ_coeff.T)
    nuclear_charges = _read_nuclear_charges(data.filename, data.mol.Natoms)
    return aobasis, occ_coeff, P_density, nuclear_charges


class _UnsupportedShellError(Exception):
    pass


def _compute_esp_at_points(aobasis, occ_coeff, P_density, nuclear_charges, xyz_bohr, Natoms,
                            Nelec, verts_bohr, esp_method, esp_spacing, esp_pad, dist_floor,
                            verbose, tag):
    Nv = verts_bohr.shape[0]

    V_nuc = np.zeros(Nv)
    for A in range(Natoms):
        dd = np.linalg.norm(verts_bohr - xyz_bohr[A, :], axis=1)
        V_nuc += nuclear_charges[A] / np.maximum(dd, 1e-6)

    used_analytic = False
    if esp_method != "numeric":
        try:
            if verbose:
                print(f"g16_draw_esp_surface ({tag}): evaluating electronic ESP term analytically (McMurchie-Davidson)...")
            t0 = time.time()
            V_elec = eval_esp_analytic(aobasis, P_density, verts_bohr, verbose)
            if verbose:
                print(f"  Analytic electronic-term evaluation: {time.time()-t0:.1f} s")
            used_analytic = True
        except ValueError as e:
            if esp_method == "analytic":
                raise
            if verbose:
                print(f"g16_draw_esp_surface ({tag}): basis has shells unsupported by the analytic ESP "
                      f"method ({e}) -- falling back to the numerical Coulomb-grid-sum method.")

    if not used_analytic:
        lo_e = xyz_bohr.min(axis=0) - esp_pad
        hi_e = xyz_bohr.max(axis=0) + esp_pad
        gxe = np.arange(lo_e[0], hi_e[0] + esp_spacing / 2, esp_spacing)
        gye = np.arange(lo_e[1], hi_e[1] + esp_spacing / 2, esp_spacing)
        gze = np.arange(lo_e[2], hi_e[2] + esp_spacing / 2, esp_spacing)
        Xe, Ye, Ze = np.meshgrid(gxe, gye, gze)
        esp_grid_pts = np.column_stack([Xe.ravel(), Ye.ravel(), Ze.ravel()])
        dVe = esp_spacing ** 3

        if verbose:
            print(f"g16_draw_esp_surface ({tag}): Coulomb-sum grid {len(gxe)} x {len(gye)} x {len(gze)} "
                  f"= {esp_grid_pts.shape[0]} points (spacing={esp_spacing:.2f} Bohr, padding={esp_pad:.1f} Bohr)")

        t0 = time.time()
        rho_e = eval_density_on_grid(aobasis, occ_coeff, esp_grid_pts)
        if verbose:
            print(f"  Coulomb-grid density evaluation: {time.time()-t0:.1f} s")

        captured_electrons = np.sum(rho_e) * dVe
        scale_factor = Nelec / captured_electrons
        if verbose:
            print(f"  Electrons captured on Coulomb grid: {captured_electrons:.2f} / {Nelec} -- rescaling density by {scale_factor:.4f}")
        rho_e = rho_e * scale_factor

        V_elec = np.zeros(Nv)
        batch_size = 100
        n_batches = (Nv + batch_size - 1) // batch_size
        t0 = time.time()
        for b in range(n_batches):
            i0, i1 = b * batch_size, min((b + 1) * batch_size, Nv)
            vb = verts_bohr[i0:i1, :]
            dx = esp_grid_pts[:, 0:1] - vb[:, 0]
            dy = esp_grid_pts[:, 1:2] - vb[:, 1]
            dz = esp_grid_pts[:, 2:3] - vb[:, 2]
            dist = np.sqrt(dx**2 + dy**2 + dz**2)
            dist = np.maximum(dist, dist_floor)
            V_elec[i0:i1] = np.sum((rho_e[:, None] * dVe) / dist, axis=0)
            if verbose and ((b + 1) % max(1, n_batches // 5) == 0 or b == n_batches - 1):
                print(f"  Coulomb sum: batch {b+1}/{n_batches} ({time.time()-t0:.1f} s elapsed)")

    return V_nuc - V_elec


def g16_draw_esp_surface(data, isovalue=0.001, grid_spacing=0.15, padding=4.0, max_vertices=3000,
                          compare_to=None, esp_method="auto", esp_grid_spacing=0.20, esp_padding=6.0,
                          distance_floor=None, save_cube="", cube_spacing=0.30, cube_padding=4.0,
                          pos_color=(0.10, 0.40, 0.85), neg_color=(0.85, 0.15, 0.10), face_alpha=None,
                          surface_style="solid", show_molecule=False, show_labels=False, atom_scale=0.35,
                          title="", ax=None, verbose=True, padding_explicit=False):
    """Renders the electrostatic potential mapped onto the density
    isosurface. See G_draw_esp_surface.m's docstring for the full
    method/accuracy notes -- names here are the same options in
    snake_case.

    Returns a dict h with key 'surf' (the coloured surface
    Poly3DCollection).
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    isoval = abs(isovalue)
    spacing = grid_spacing
    pad = padding
    compare_mode = compare_to is not None
    esp_method = esp_method.lower()
    if distance_floor is None:
        distance_floor = 0.5 * esp_grid_spacing
    face_alpha_explicit = face_alpha is not None
    if face_alpha is None:
        face_alpha = 1.0 if surface_style == "solid" else 0.55
    elif surface_style == "transparent" and not face_alpha_explicit:
        face_alpha = 0.55

    if show_molecule and surface_style != "grid":
        if not face_alpha_explicit and surface_style == "solid":
            face_alpha = 0.6
            if verbose:
                print(f"g16_draw_esp_surface: show_molecule=True, so face_alpha defaults to {face_alpha:.2g} instead of opaque 1.0.")
        elif face_alpha_explicit and face_alpha > 0.9:
            print(f"g16_draw_esp_surface: WARNING show_molecule=True but face_alpha={face_alpha:.2g} is nearly opaque.")

    aobasis, occ_coeff, P_density, nuclear_charges = _validate_and_prepare(data, "data")
    if compare_mode:
        aobasis2, occ_coeff2, P_density2, nuclear_charges2 = _validate_and_prepare(compare_to, "compare_to")

    xyz_bohr = np.asarray(data.xyz_bohr)

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
                raise ValueError(f"g16_draw_esp_surface: requested surface-shape grid has {Ngrid} points, too large.")
            print(f"g16_draw_esp_surface: stopped auto-expanding padding at {pad:.1f} Bohr.")
            break

        if verbose:
            print(f"g16_draw_esp_surface: surface-shape grid {len(gx)} x {len(gy)} x {len(gz)} "
                  f"= {Ngrid} points (spacing={spacing:.2f} Bohr, padding={pad:.1f} Bohr)")

        X, Y, Z = np.meshgrid(gx, gy, gz)
        grid_pts_bohr = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

        t0 = time.time()
        V = eval_density_on_grid(aobasis, occ_coeff, grid_pts_bohr).reshape(X.shape)
        if verbose:
            print(f"  Density evaluation: {time.time()-t0:.1f} s")

        Xa, Ya, Za = gx * _A0, gy * _A0, gz * _A0
        verts, faces, normals = marching_cubes_ang(Xa, Ya, Za, V, isoval)
        clipped = is_clipped(verts, Xa, Ya, Za, 1.5 * spacing * _A0)

        if not clipped or padding_explicit or attempt >= max_attempts:
            break
        if verbose:
            print(f"  Isosurface touches the grid boundary -- expanding padding {pad:.1f} -> {pad*pad_growth:.1f} Bohr.")
        pad *= pad_growth

    if clipped:
        print(f"g16_draw_esp_surface: WARNING density isosurface appears clipped (padding={pad:.1f} Bohr).")
    if verts is None:
        raise ValueError(f"g16_draw_esp_surface: no density isosurface found at isovalue={isoval:.4g}.")

    n_full = verts.shape[0]
    dec = decimate_mesh(verts, faces, max_vertices, normals=normals)
    verts, faces, normals = dec["vertices"], dec["faces"], dec["normals"]
    if verbose and verts.shape[0] < n_full:
        print(f"  Mesh decimated: {n_full} -> {verts.shape[0]} vertices")

    verts_bohr = verts / _A0

    V_esp1 = _compute_esp_at_points(aobasis, occ_coeff, P_density, nuclear_charges, xyz_bohr,
                                     data.mol.Natoms, data.Nelec, verts_bohr, esp_method,
                                     esp_grid_spacing, esp_padding, distance_floor, verbose, "data")
    if compare_mode:
        V_esp2 = _compute_esp_at_points(aobasis2, occ_coeff2, P_density2, nuclear_charges2,
                                         np.asarray(compare_to.xyz_bohr), compare_to.mol.Natoms,
                                         compare_to.Nelec, verts_bohr, esp_method, esp_grid_spacing,
                                         esp_padding, distance_floor, verbose, "compare_to")
        V_esp = V_esp2 - V_esp1
    else:
        V_esp = V_esp1

    if save_cube:
        max_cube_grid = 2_000_000
        lo_c = xyz_bohr.min(axis=0) - cube_padding
        hi_c = xyz_bohr.max(axis=0) + cube_padding
        gxc = np.arange(lo_c[0], hi_c[0] + cube_spacing / 2, cube_spacing)
        gyc = np.arange(lo_c[1], hi_c[1] + cube_spacing / 2, cube_spacing)
        gzc = np.arange(lo_c[2], hi_c[2] + cube_spacing / 2, cube_spacing)
        Ngrid_cube = len(gxc) * len(gyc) * len(gzc)
        if Ngrid_cube > max_cube_grid:
            raise ValueError(f"g16_draw_esp_surface: save_cube grid has {Ngrid_cube} points, too large (cap {max_cube_grid}).")
        if verbose:
            print(f"g16_draw_esp_surface: save_cube grid {len(gxc)} x {len(gyc)} x {len(gzc)} = {Ngrid_cube} points.")
        Xc, Yc, Zc = np.meshgrid(gxc, gyc, gzc)
        cube_pts_bohr = np.column_stack([Xc.ravel(), Yc.ravel(), Zc.ravel()])
        Vc1 = _compute_esp_at_points(aobasis, occ_coeff, P_density, nuclear_charges, xyz_bohr,
                                      data.mol.Natoms, data.Nelec, cube_pts_bohr, esp_method,
                                      esp_grid_spacing, esp_padding, distance_floor, verbose, "data, cube grid")
        if compare_mode:
            Vc2 = _compute_esp_at_points(aobasis2, occ_coeff2, P_density2, nuclear_charges2,
                                          np.asarray(compare_to.xyz_bohr), compare_to.mol.Natoms,
                                          compare_to.Nelec, cube_pts_bohr, esp_method, esp_grid_spacing,
                                          esp_padding, distance_floor, verbose, "compare_to, cube grid")
            Vcube = (Vc2 - Vc1).reshape(Xc.shape)
            title_line = "ESP difference (compare_to - data) -- G_Utility export"
        else:
            Vcube = Vc1.reshape(Xc.shape)
            title_line = "Electrostatic potential -- G_Utility export"
        comment_line = f"V(r), Hartree/e, grid spacing={cube_spacing:.3f} Bohr"
        write_cube_file(save_cube, title_line, comment_line, data.mol.Z, xyz_bohr, gxc, gyc, gzc, Vcube)
        if verbose:
            print(f"  Cube file written: {save_cube}")

    if ax is None:
        fig = plt.figure()
        fig.canvas.manager.set_window_title("ESP surface")
        ax = fig.add_subplot(111, projection="3d")

    cmap = diverging_cmap(neg_color, pos_color)
    vmax = pctile_local(np.abs(V_esp), 99)
    if vmax == 0:
        vmax = 1.0

    face_vals = face_average_values(V_esp, faces)
    face_colors = cmap((face_vals + vmax) / (2 * vmax))

    if surface_style == "grid":
        h_surf = Poly3DCollection(verts[faces], facecolors=(0, 0, 0, 0),
                                   edgecolors=face_colors, linewidths=0.6, alpha=1.0)
    else:
        h_surf = Poly3DCollection(verts[faces], facecolors=face_colors,
                                   edgecolors=(0, 0, 0, 0), alpha=face_alpha, shade=False)
    ax.add_collection3d(h_surf)

    import matplotlib.cm as mcm
    import matplotlib.colors as mcolors
    sm = mcm.ScalarMappable(norm=mcolors.Normalize(vmin=-vmax, vmax=vmax), cmap=cmap)
    cb = plt.colorbar(sm, ax=ax)
    cb.set_label("ESP difference, compare_to - data (Hartree/e)" if compare_mode else "Electrostatic potential (Hartree/e)")

    finish_3d_axes(ax, xyz_bohr * _A0, verts)

    if show_molecule:
        draw_molecule_or_fallback(ax, data.mol, show_labels, atom_scale, "g16_draw_esp_surface")

    if not title:
        title = (f"ESP difference on density surface (iso={isoval:.3g})" if compare_mode
                  else f"ESP on density surface (iso={isoval:.3g})")
    ax.set_title(title, fontsize=11)

    return {"surf": h_surf}
