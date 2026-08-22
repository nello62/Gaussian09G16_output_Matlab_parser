import os
import subprocess
import tempfile

import numpy as np

from ._common import Struct
from .draw_molecule import _classify_bond_order, g16_draw_molecule
from .get_bond_length import g16_get_bond_length


def _content_bbox(img):
    """Tight (row_min, row_max, col_min, col_max) bounding box (inclusive,
    0-based) of the non-background pixels in a rendered RGB(A) frame
    (img: HxWx3 or HxWx4 uint8). The figure background is plain white, so
    any pixel that is not near-white is molecule/label/title content.
    Falls back to the full image if nothing is found (e.g. an all-white
    frame), so cropping degrades to a no-op rather than erroring. Same
    algorithm as G16_animate_mode.m's local_content_bbox.
    """
    is_background = np.all(img[..., :3] >= 250, axis=-1)
    rows, cols = np.where(~is_background)
    if rows.size == 0:
        return (0, img.shape[0] - 1, 0, img.shape[1] - 1)
    return (int(rows.min()), int(rows.max()), int(cols.min()), int(cols.max()))


def _union_bbox(b1, b2):
    return (min(b1[0], b2[0]), max(b1[1], b2[1]), min(b1[2], b2[2]), max(b1[3], b2[3]))


def _pad_bbox(bbox, margin, img_shape):
    """Expands bbox by margin pixels on every side, clamped to img_shape
    (H, W, ...), and nudges width/height to be even -- required by
    libx264/yuv420p, which rejects odd frame dimensions. Same algorithm as
    G16_animate_mode.m's local_pad_bbox.
    """
    row_min, row_max, col_min, col_max = bbox
    h, w = img_shape[0], img_shape[1]
    row_min = max(0, row_min - margin)
    row_max = min(h - 1, row_max + margin)
    col_min = max(0, col_min - margin)
    col_max = min(w - 1, col_max + margin)

    if (row_max - row_min + 1) % 2 != 0:
        if row_max < h - 1:
            row_max += 1
        else:
            row_min -= 1
    if (col_max - col_min + 1) % 2 != 0:
        if col_max < w - 1:
            col_max += 1
        else:
            col_min -= 1
    return (row_min, row_max, col_min, col_max)


def g16_animate_mode(mol, nm, mode_idx, filename=None, scale=1.5, flip_sign=False,
                      atom_scale=0.35, bond_tol=1.30, show_labels=False,
                      frames_per_cycle=30, n_cycles=2, fps=20, view=None,
                      dpi=150, progress_callback=None,
                      crop=True, crop_margin=12, show_title=False):
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
    dpi : int — figure resolution for the saved video (default 150; the
        matplotlib figure default of 100 combined with the default figure
        size gives only 640x480, encoded at a low bitrate that looks
        noticeably compressed/blocky -- 150 gives a sharper 960x720 video)
    view : tuple (azim, elev) in degrees, optional — starting camera
        orientation (default None = matplotlib's default 3D view). Pass
        (ax.azim, ax.elev) from a figure you have already rotated
        interactively.
    progress_callback : callable(current_frame, total_frames), optional —
        forwarded to matplotlib's FuncAnimation.save(), called after each
        frame is written. Rendering + encoding a full animation can take
        well over a minute for a large molecule with no other feedback
        otherwise; useful for a caller (e.g. a GUI) to show progress
        instead of appearing frozen.
    crop : bool — tightly crop the video around the molecule instead of
        the whole (mostly blank) figure canvas (default True). The crop
        rectangle is computed once, from the two oscillation extremes
        (phase +1 and -1), then held fixed for every frame, so the video
        does not jitter or clip atoms at any phase. Implemented as an
        ffmpeg `crop` filter pass over the full-resolution video (ffmpeg
        is already required for MP4 output), rather than per-frame pixel
        manipulation in Python.
    crop_margin : int — blank margin kept around the molecule when crop is
        True, in pixels (default 12)
    show_title : bool — draw the filename/mode/frequency title text in
        every frame (default False). Left off by default so the title's
        own bounding box does not force a taller crop than the molecule
        itself needs.

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

    # Force equal aspect ratio (same cube-shaped box on all three axes),
    # same as g16_draw_molecule's own _set_axes_equal. Without this, a
    # molecule whose bounding box is not roughly cubic (e.g. a mostly
    # planar/aromatic one, common here) gets one axis visibly compressed
    # relative to the others -- CPK spheres render as flattened ellipses,
    # not obvious from every camera angle but very visible once rotated
    # to look more edge-on along the compressed axis.
    _mid = [(lo + hi) / 2 for lo, hi in (xlim, ylim, zlim)]
    _radius = max(hi - lo for lo, hi in (xlim, ylim, zlim)) / 2
    xlim = (_mid[0] - _radius, _mid[0] + _radius)
    ylim = (_mid[1] - _radius, _mid[1] + _radius)
    zlim = (_mid[2] - _radius, _mid[2] + _radius)

    # Fixed bond list (and bond order) from the equilibrium geometry, so
    # bonds do not appear/disappear or flicker between single/double/
    # triple frame to frame as instantaneous distances oscillate across
    # the bond_tol threshold (g16_draw_molecule's default distance-based
    # detection would otherwise re-evaluate both on every frame).
    bond_table = g16_get_bond_length(mol, tolerance=bond_tol, include_h=True)
    orders = [
        _classify_bond_order(s1, s2, d)
        for s1, s2, d in zip(bond_table["Sym1"], bond_table["Sym2"], bond_table["Distance_Ang"])
    ]
    bond_list = bond_table[["Atom1", "Atom2"]].to_numpy() - 1
    bond_list = np.column_stack([bond_list, orders])

    freq_str = f"Mode {mode_idx} - {nm.freq[i0]:.1f} cm$^{{-1}}$"
    if show_title:
        fname = os.path.splitext(os.path.basename(mol.filename))[0] if getattr(mol, "filename", None) else ""
        title = f"{fname}\n{freq_str}" if fname else freq_str
    else:
        # A single space renders no visible ink, so it is excluded from
        # the content bounding box used for cropping, unlike an empty
        # string (which g16_draw_molecule replaces with a default title).
        title = " "

    fig = plt.figure()
    fig.set_dpi(dpi)  # match the dpi ani.save() will use, so a probe
                       # frame's pixel bbox applies unchanged to the video
    ax = fig.add_subplot(111, projection="3d")

    total_frames = frames_per_cycle * n_cycles
    mol_frame = Struct(**vars(mol))

    def render_phase(phase):
        ax.clear()
        mol_frame.xyz = mol.xyz + phase * U_scaled
        g16_draw_molecule(mol_frame, ax=ax, atom_scale=atom_scale, bond_tol=bond_tol,
                           show_labels=show_labels, show_legend=False, title=title,
                           bond_list=bond_list)
        if view is not None:
            ax.view_init(elev=view[1], azim=view[0])
        ax.set_xlim3d(xlim)
        ax.set_ylim3d(ylim)
        ax.set_zlim3d(zlim)

    def update(k):
        render_phase(np.sin(2 * np.pi * k / frames_per_cycle))
        return []

    # -----------------------------------------------------------------
    # Crop rectangle: fixed axis limits mean the plotted box occupies the
    # same pixel rectangle in every frame regardless of oscillation phase,
    # so it is enough to measure it once here (from the two oscillation
    # extremes, which bound the whole animation) and reuse it for every
    # frame -- this is what removes the large blank margin around the
    # molecule without ever clipping an atom at any phase.
    # -----------------------------------------------------------------
    crop_box = None
    if crop:
        render_phase(1)
        fig.canvas.draw()
        buf_plus = np.asarray(fig.canvas.buffer_rgba())
        render_phase(-1)
        fig.canvas.draw()
        buf_minus = np.asarray(fig.canvas.buffer_rgba())
        bbox = _union_bbox(_content_bbox(buf_plus), _content_bbox(buf_minus))
        crop_box = _pad_bbox(bbox, crop_margin, buf_plus.shape)

    # bitrate=-1 (do not force a fixed bitrate) + a CRF-based quality
    # target via extra_args gives a much sharper result than
    # FFMpegWriter's own low-bitrate default (previously ~160 kbps at
    # 640x480, visibly blocky); yuv420p keeps the output playable in
    # QuickTime and most other players, which can choke on matplotlib's
    # raw output pixel format otherwise.
    writer = animation.FFMpegWriter(
        fps=fps, bitrate=-1,
        extra_args=["-crf", "18", "-preset", "slow", "-pix_fmt", "yuv420p"],
    )
    ani = animation.FuncAnimation(fig, update, frames=total_frames, blit=False)

    if crop_box is None:
        ani.save(filename, writer=writer, dpi=dpi, progress_callback=progress_callback)
    else:
        # Render the full, uncropped video to a temp file, then crop it
        # with a single ffmpeg pass -- simpler and far more robust than
        # cropping matplotlib's raw per-frame pixel buffers ourselves,
        # and ffmpeg is already a hard requirement for MP4 output at all.
        row_min, row_max, col_min, col_max = crop_box
        crop_w = col_max - col_min + 1
        crop_h = row_max - row_min + 1
        tmp_fd, tmp_path = tempfile.mkstemp(suffix=".mp4")
        os.close(tmp_fd)
        try:
            ani.save(tmp_path, writer=writer, dpi=dpi, progress_callback=progress_callback)
            result = subprocess.run(
                ["ffmpeg", "-y", "-i", tmp_path,
                 "-vf", f"crop={crop_w}:{crop_h}:{col_min}:{row_min}",
                 "-c:v", "libx264", "-crf", "18", "-preset", "slow", "-pix_fmt", "yuv420p",
                 filename],
                capture_output=True, text=True,
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"g16_animate_mode: ffmpeg crop pass failed (exit {result.returncode}):\n{result.stderr}"
                )
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    plt.close(fig)

    print(f"g16_animate_mode: animation saved to {filename} ({total_frames} frames, {fps} fps)")
    return filename
