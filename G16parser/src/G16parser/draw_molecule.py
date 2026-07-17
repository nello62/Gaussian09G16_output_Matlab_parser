import os

import numpy as np

from ._common import get_color, get_radius, DEFAULT_COLOR


def _set_axes_equal(ax):
    """Equal 3D aspect ratio workaround (matplotlib < 3.3 has no set_box_aspect)."""
    limits = np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
    middle = limits.mean(axis=1)
    radius = 0.5 * max(abs(limits[:, 1] - limits[:, 0]))
    ax.set_xlim3d([middle[0] - radius, middle[0] + radius])
    ax.set_ylim3d([middle[1] - radius, middle[1] + radius])
    ax.set_zlim3d([middle[2] - radius, middle[2] + radius])


def _draw_cartesian_axes(ax, origin, axes_length):
    colors = [(0.85, 0.10, 0.10), (0.10, 0.65, 0.10), (0.10, 0.10, 0.85)]
    labels = ["X", "Y", "Z"]
    dirs = np.eye(3) * axes_length
    for a in range(3):
        d = dirs[a]
        ax.quiver(origin[0], origin[1], origin[2], d[0], d[1], d[2],
                  color=colors[a], linewidth=1.5)
        tip = origin + d
        ax.text(tip[0], tip[1], tip[2], f"  {labels[a]}", color=colors[a],
                fontsize=10, fontweight="bold")


def g16_draw_molecule(mol, atom_scale=0.35, bond_tol=1.30, show_labels=True,
                       show_legend=True, title="", bg_color=(0.95, 0.95, 0.95),
                       ax=None, show_axes=False, axes_length=None, bond_list=None):
    """Renders a 3D CPK ball-and-stick model from the mol struct (matplotlib,
    static — no interactive rotation widget, unlike the MATLAB original's
    ``rotate3d``; use the normal matplotlib 3D toolbar/mouse controls
    instead when shown in an interactive backend).

    Parameters mirror G16_draw_molecule.m — see its docstring for details.
    bond_list : array-like (Nbonds, 2) of 0-based atom-index pairs,
        optional — draws exactly these bonds instead of auto-detecting
        from bond_tol. Useful to keep a fixed bond topology across a
        series of frames where atoms move (e.g. g16_animate_mode), so
        bonds do not appear/disappear as instantaneous distances cross
        the bond_tol threshold.

    Returns
    -------
    ax : the matplotlib 3D axes used
    """
    import matplotlib.pyplot as plt  # noqa: F401 (ensures a backend is selected)

    if not hasattr(mol, "symbols") or not hasattr(mol, "xyz"):
        raise ValueError("g16_draw_molecule: mol must be the struct returned by g16_structure.")

    if not title:
        if getattr(mol, "filename", None):
            fn = os.path.splitext(os.path.basename(mol.filename))[0]
            title = fn
        else:
            title = "Molecule"

    if ax is None:
        fig = plt.figure()
        fig.canvas.manager.set_window_title(title)
        ax = fig.add_subplot(111, projection="3d")

    ax.set_facecolor(bg_color)

    # -----------------------------------------------------------------
    # Bonds
    # -----------------------------------------------------------------
    xyz = mol.xyz
    symbols = mol.symbols
    if bond_list is None:
        for i in range(mol.Natoms):
            ri = get_radius(symbols[i])
            for j in range(i + 1, mol.Natoms):
                rj = get_radius(symbols[j])
                d = np.linalg.norm(xyz[i] - xyz[j])
                if d < (ri + rj) * bond_tol:
                    p1, p2 = xyz[i], xyz[j]
                    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                            color=(0.5, 0.5, 0.5), linewidth=2.0)
    else:
        for i, j in bond_list:
            p1, p2 = xyz[i], xyz[j]
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                    color=(0.5, 0.5, 0.5), linewidth=2.0)

    # -----------------------------------------------------------------
    # Atoms (spheres), heavy elements first, H last (for the legend order)
    # -----------------------------------------------------------------
    seen = []
    for s in symbols:
        if s not in seen:
            seen.append(s)
    heavy = [s for s in seen if s != "H"]
    syms_ordered = heavy + (["H"] if "H" in seen else [])

    u, v = np.mgrid[0:2 * np.pi:24j, 0:np.pi:16j]
    xs, ys, zs = np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v)

    legend_handles, legend_labels = [], []

    for sym in syms_ordered:
        clr = get_color(sym)
        r = get_radius(sym) * atom_scale
        idx_atoms = [i for i, s in enumerate(symbols) if s == sym]

        for i in idx_atoms:
            cx, cy, cz = xyz[i]
            ax.plot_surface(xs * r + cx, ys * r + cy, zs * r + cz,
                             color=clr, shade=True, linewidth=0, antialiased=False)
            if show_labels:
                off = r * 1.4
                dim_clr = tuple(c * 0.7 for c in clr)
                ax.text(cx + off, cy + off, cz + off, f"{sym}{i+1}", fontsize=7,
                        color=dim_clr, ha="left")

        proxy, = ax.plot([], [], [], "o", color=clr, markersize=8)
        legend_handles.append(proxy)
        legend_labels.append(sym)

    if show_legend:
        ax.legend(legend_handles, legend_labels, loc="upper left", fontsize=9,
                  frameon=False)

    ax.set_title(title, fontsize=11)
    ax.set_axis_off()
    _set_axes_equal(ax)

    if show_axes:
        xlim, ylim, zlim = ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()
        corner = np.array([xlim[0], ylim[0], zlim[0]])
        diag_len = np.linalg.norm([xlim[1] - xlim[0], ylim[1] - ylim[0], zlim[1] - zlim[0]])
        length = axes_length if axes_length else (diag_len * 0.20 or 2.0)
        _draw_cartesian_axes(ax, corner, length)

    try:
        ax.set_proj_type("persp")
    except Exception:
        pass

    return ax
