"""2D cutting-plane helpers for 'mode'='contour': plane basis construction
(PCA best-fit or coordinate planes) and a lightweight atom/bond sketch.
Direct port of g_plane_basis.m / g_draw_2d_atoms.m.
"""

import numpy as np

_CPK_COLORS = {
    "H": (0.60, 0.80, 1.00), "C": (0.30, 0.30, 0.30), "N": (0.10, 0.30, 0.90),
    "O": (0.90, 0.10, 0.10), "F": (0.20, 0.80, 0.20), "P": (1.00, 0.50, 0.00),
    "S": (1.00, 0.85, 0.00), "Cl": (0.20, 0.85, 0.20), "Br": (0.55, 0.20, 0.10),
    "I": (0.45, 0.00, 0.65), "Au": (1.00, 0.82, 0.14),
}
_COV_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07,
    "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39, "Au": 1.36,
}
_DEFAULT_COLOR = (0.65, 0.20, 0.80)
_DEFAULT_RADIUS = 0.80
_BOND_COLOR = (0.45, 0.45, 0.45)


def plane_basis(xyz_bohr, plane_choice):
    """Returns center [3], orthonormal in-plane basis u, v [3], and plane
    normal n [3] for 'mode'='contour'. plane_choice: 'auto' (PCA best-fit
    plane through all atoms via SVD) | 'xy' | 'xz' | 'yz'.
    """
    xyz_bohr = np.asarray(xyz_bohr)
    center = xyz_bohr.mean(axis=0)
    if plane_choice == "auto":
        if xyz_bohr.shape[0] < 3:
            u, v, n = np.array([1., 0, 0]), np.array([0., 1, 0]), np.array([0., 0, 1])
        else:
            _, _, Vt = np.linalg.svd(xyz_bohr - center, full_matrices=False)
            u, v, n = Vt[0], Vt[1], Vt[2]
    elif plane_choice == "xy":
        u, v, n = np.array([1., 0, 0]), np.array([0., 1, 0]), np.array([0., 0, 1])
    elif plane_choice == "xz":
        u, v, n = np.array([1., 0, 0]), np.array([0., 0, 1]), np.array([0., 1, 0])
    elif plane_choice == "yz":
        u, v, n = np.array([0., 1, 0]), np.array([0., 0, 1]), np.array([1., 0, 0])
    else:
        raise ValueError(f"plane_basis: plane_choice must be 'auto'/'xy'/'xz'/'yz', got {plane_choice!r}")
    return center, u, v, n


def draw_2d_atoms(ax, symbols, atom_xy, xyz_ang):
    """Lightweight 2D sketch of atoms (CPK-coloured markers + element/
    index labels) and simple distance-based bonds, projected onto a
    contour cutting plane. atom_xy [Nat,2] already-projected in-plane
    positions (Angstrom); xyz_ang [Nat,3] true 3D positions (Angstrom),
    used only for the bond-length distance test.
    """
    Nat = len(symbols)

    def radius(sym):
        return _COV_RADII.get(sym, _DEFAULT_RADIUS)

    def color(sym):
        return _CPK_COLORS.get(sym, _DEFAULT_COLOR)

    for i in range(Nat):
        ri = radius(symbols[i])
        for j in range(i + 1, Nat):
            rj = radius(symbols[j])
            d = np.linalg.norm(xyz_ang[i] - xyz_ang[j])
            if d < (ri + rj) * 1.30:
                ax.plot([atom_xy[i, 0], atom_xy[j, 0]], [atom_xy[i, 1], atom_xy[j, 1]],
                        "-", color=_BOND_COLOR, linewidth=1.5)

    for i in range(Nat):
        clr = color(symbols[i])
        ax.plot(atom_xy[i, 0], atom_xy[i, 1], "o", markerfacecolor=clr,
                markeredgecolor="k", markersize=7)
        dark = tuple(c * 0.6 for c in clr)
        ax.text(atom_xy[i, 0], atom_xy[i, 1], f"  {symbols[i]}{i+1}", fontsize=7, color=dark)
