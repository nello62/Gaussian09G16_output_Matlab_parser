import os
import warnings

import numpy as np

from .draw_molecule import g16_draw_molecule


def g16_draw_mode(mol, nm, mode_idx, scale=1.5, arrow_color=(1.0, 0.4, 0.1),
                  atom_scale=0.35, bond_tol=1.30, show_labels=False, flip_sign=False):
    """Visualises a normal mode on a 3D molecular structure (matplotlib
    quiver arrows in place of MATLAB's hand-built cone-tipped arrows —
    same information, simpler static rendering).

    Parameters
    ----------
    mol : Struct — from g16_structure
    nm : Struct — from g16_nmodes
    mode_idx : int — 1-based mode index into nm.freq
    scale : float — arrow length scale (default 1.5)
    arrow_color : tuple — default (1.0, 0.4, 0.1)
    atom_scale, bond_tol, show_labels : see g16_draw_molecule
    flip_sign : bool — invert all arrow directions (default False); normal
        mode eigenvectors have an arbitrary overall sign, so a mode and its
        180-degree-phase-shifted twin are physically identical.

    Returns
    -------
    ax : the matplotlib 3D axes used
    """
    import matplotlib.pyplot as plt

    if mode_idx < 1 or mode_idx > nm.Nmodes:
        raise ValueError(f"g16_draw_mode: mode index {mode_idx} is out of range [1, {nm.Nmodes}]")
    if mol.Natoms != nm.Natoms:
        raise ValueError(f"g16_draw_mode: mol.Natoms ({mol.Natoms}) does not match nm.Natoms ({nm.Natoms})")

    i0 = mode_idx - 1
    freq_str = f"Mode {mode_idx}  -  {nm.freq[i0]:.1f} cm-1"
    if nm.has_Raman:
        freq_str += f"   IR={nm.IR[i0]:.1f}   Raman={nm.Raman[i0]:.1f}"

    fig = plt.figure()
    fig.canvas.manager.set_window_title(f"Mode {mode_idx} - {nm.freq[i0]:.1f} cm-1")
    ax = fig.add_subplot(111, projection="3d")

    g16_draw_molecule(mol, ax=ax, atom_scale=atom_scale, bond_tol=bond_tol,
                       show_labels=show_labels, show_legend=False, title=freq_str)

    U = nm.disp[:, :, i0].copy()
    if flip_sign:
        U = -U

    norms_i = np.linalg.norm(U, axis=1)
    max_norm = norms_i.max()
    if max_norm == 0:
        warnings.warn(f"g16_draw_mode: zero displacement vectors for mode {mode_idx}")
        return ax

    U_scaled = U / max_norm * scale
    thresh = 0.05

    for i in range(mol.Natoms):
        if norms_i[i] / max_norm < thresh:
            continue
        x0, y0, z0 = mol.xyz[i]
        dx, dy, dz = U_scaled[i]
        ax.quiver(x0, y0, z0, dx, dy, dz, color=arrow_color, linewidth=2.0,
                  arrow_length_ratio=0.3)

    fname = os.path.splitext(os.path.basename(mol.filename))[0] if getattr(mol, "filename", None) else ""
    ax.set_title(f"{fname}\n{freq_str}" if fname else freq_str, fontsize=10)

    return ax
