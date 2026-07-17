import os

import numpy as np

from ._common import Struct
from .draw_molecule import g16_draw_molecule
from .get_bond_length import g16_get_bond_length


def g16_animate_mode(mol, nm, mode_idx, filename=None, scale=1.5, flip_sign=False,
                      atom_scale=0.35, bond_tol=1.30, show_labels=False,
                      frames_per_cycle=30, n_cycles=2, fps=20, view=None):
    """Exports an MP4 animation of a vibrational mode.

    Python port of G16_animate_mode.m: oscillates the molecule along the
    mode's displacement vector (equilibrium +/- amplitude, like GaussView's
    mode animations) and saves the result as an MP4 video via matplotlib's
    FFMpegWriter.

    Requires ffmpeg to be installed and on PATH — matplotlib does not
    bundle a video encoder itself:
        macOS (Homebrew):  brew install ffmpeg
        Ubuntu/Debian:      sudo apt install ffmpeg

    Parameters
    ----------
    mol : Struct — from g16_structure
    nm : Struct — from g16_nmodes
    mode_idx : int — 1-based mode index into nm.freq
    filename : str, optional — output path (default: "<source>_mode<N>.mp4";
        ".mp4" is appended if missing)
    scale : float — displacement amplitude scale, same meaning as
        g16_draw_mode's scale (default 1.5)
    flip_sign : bool — invert the displacement direction (default False)
    atom_scale, bond_tol, show_labels : see g16_draw_molecule
    frames_per_cycle : int — frames per oscillation period (default 30)
    n_cycles : int — number of periods rendered (default 2)
    fps : int — video frame rate (default 20)
    view : tuple (azim, elev) in degrees, optional — starting camera
        orientation (default None = matplotlib's default 3D view). Pass
        (ax.azim, ax.elev) from a figure you have already rotated
        interactively.

    Returns
    -------
    filename : str — the path written to.
    """
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt

    if mode_idx < 1 or mode_idx > nm.Nmodes:
        raise ValueError(f"g16_animate_mode: mode index {mode_idx} is out of range [1, {nm.Nmodes}]")
    if mol.Natoms != nm.Natoms:
        raise ValueError(f"g16_animate_mode: mol.Natoms ({mol.Natoms}) does not match nm.Natoms ({nm.Natoms})")

    if not filename:
        fn = os.path.splitext(os.path.basename(mol.filename))[0] if getattr(mol, "filename", None) else "molecule"
        filename = f"{fn}_mode{mode_idx}.mp4"
    if not filename.lower().endswith(".mp4"):
        filename += ".mp4"

    i0 = mode_idx - 1
    U = nm.disp[:, :, i0].copy()
    if flip_sign:
        U = -U
    norms_i = np.linalg.norm(U, axis=1)
    max_norm = norms_i.max()
    if max_norm == 0:
        raise ValueError(f"g16_animate_mode: zero displacement vectors for mode {mode_idx}")
    U_scaled = U / max_norm * scale

    # Fixed axis limits across the whole oscillation, so the camera/box
    # does not jitter frame to frame.
    pad = 1.0
    extreme = np.vstack([mol.xyz - np.abs(U_scaled), mol.xyz + np.abs(U_scaled)])
    xlim = (extreme[:, 0].min() - pad, extreme[:, 0].max() + pad)
    ylim = (extreme[:, 1].min() - pad, extreme[:, 1].max() + pad)
    zlim = (extreme[:, 2].min() - pad, extreme[:, 2].max() + pad)

    # Fixed bond list from the equilibrium geometry, so bonds do not
    # appear/disappear frame to frame as instantaneous distances oscillate
    # across the bond_tol threshold (g16_draw_molecule's default
    # distance-based detection would otherwise re-evaluate connectivity
    # on every frame).
    bond_table = g16_get_bond_length(mol, tolerance=bond_tol, include_h=True)
    bond_list = bond_table[["Atom1", "Atom2"]].to_numpy() - 1

    freq_str = f"Mode {mode_idx} - {nm.freq[i0]:.1f} cm$^{{-1}}$"
    fname = os.path.splitext(os.path.basename(mol.filename))[0] if getattr(mol, "filename", None) else ""
    title = f"{fname}\n{freq_str}" if fname else freq_str

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    total_frames = frames_per_cycle * n_cycles
    mol_frame = Struct(**vars(mol))

    def update(k):
        ax.clear()
        phase = np.sin(2 * np.pi * k / frames_per_cycle)
        mol_frame.xyz = mol.xyz + phase * U_scaled
        g16_draw_molecule(mol_frame, ax=ax, atom_scale=atom_scale, bond_tol=bond_tol,
                           show_labels=show_labels, show_legend=False, title=title,
                           bond_list=bond_list)
        if view is not None:
            ax.view_init(elev=view[1], azim=view[0])
        ax.set_xlim3d(xlim)
        ax.set_ylim3d(ylim)
        ax.set_zlim3d(zlim)
        return []

    ani = animation.FuncAnimation(fig, update, frames=total_frames, blit=False)
    ani.save(filename, writer=animation.FFMpegWriter(fps=fps))
    plt.close(fig)

    print(f"g16_animate_mode: animation saved to {filename} ({total_frames} frames, {fps} fps)")
    return filename
