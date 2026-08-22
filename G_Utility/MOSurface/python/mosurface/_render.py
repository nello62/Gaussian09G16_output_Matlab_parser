"""Rendering helpers shared by all three grid-evaluated
g16_draw_*_surface functions (mo/density/esp): SurfaceStyle -> matplotlib
Poly3DCollection kwargs, diverging colormap construction, 3D axes
framing, and the CPK molecule overlay dispatch.
"""

import numpy as np


def style_kwargs(style, color, face_alpha):
    """Poly3DCollection kwargs for SurfaceStyle: 'grid' draws only the
    mesh wireframe (no face fill, edges coloured); 'solid'/'transparent'
    fill faces with a flat shaded colour. matplotlib's Poly3DCollection
    has no true per-vertex Gouraud shading like MATLAB's PATCH --
    shade=True gives simple per-face directional shading instead, the
    closest native equivalent.
    """
    if style == "grid":
        return {"facecolors": (0, 0, 0, 0), "edgecolors": [color], "linewidths": 0.5, "alpha": 1.0, "shade": False}
    return {"facecolors": [color], "edgecolors": (0, 0, 0, 0), "alpha": face_alpha, "shade": True}


def colorby_style_kwargs(style, face_alpha):
    """Like style_kwargs, but for a per-vertex/per-face coloured surface
    ('ColorBy' in g16_draw_cube_surface, or the ESP map in
    g16_draw_esp_surface): colours come from facecolors set by the
    caller (per-face, averaged from per-vertex data -- matplotlib's
    Poly3DCollection has no interpolated per-vertex FaceColor the way
    MATLAB's PATCH('FaceColor','interp') does), so this only supplies
    the alpha/edge/shade behaviour, not colours.
    """
    if style == "grid":
        return {"linewidths": 0.5, "alpha": 1.0, "shade": False}
    return {"edgecolors": (0, 0, 0, 0), "alpha": face_alpha, "shade": True}


def diverging_cmap(neg_color, pos_color, n=256):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("diverging", [neg_color, (1, 1, 1), pos_color], N=n)


def sequential_cmap(pos_color, n=256):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("sequential", [(1, 1, 1), pos_color], N=n)


def finish_3d_axes(ax, xyz_ang, *vertex_arrays):
    """Sets equal-aspect, tight limits around both the isosurface(s) and
    the atoms -- the matplotlib equivalent of MATLAB's 'axis equal'/'axis
    tight' pair used by every G_draw_*_surface caller.
    """
    pts = [xyz_ang]
    for v in vertex_arrays:
        if v is not None:
            pts.append(v)
    allpts = np.vstack(pts)
    lo, hi = allpts.min(axis=0), allpts.max(axis=0)
    center = (lo + hi) / 2
    radius = max((hi - lo).max() / 2, 1e-6)
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()


def draw_molecule_or_fallback(ax, mol, show_labels, atom_scale, fname_prefix="g16_draw_surface"):
    """Overlays the CPK ball-and-stick model via G16parser's
    g16_draw_molecule -- the Python equivalent of the MATLAB original's
    G09/G16_draw_molecule dispatch (there is only ever one Python
    implementation, so no version dispatch is needed).
    """
    try:
        from G16parser import g16_draw_molecule
    except ImportError as e:
        raise ImportError(
            f"{fname_prefix}: show_molecule=True requires G16parser to be "
            "importable (pip install -e G16parser/)."
        ) from e
    g16_draw_molecule(mol, ax=ax, show_labels=show_labels, atom_scale=atom_scale,
                       show_legend=True, title="")


def face_average_values(vertex_values, faces):
    """Per-face colour value = mean of its 3 vertices' scalar values --
    the practical matplotlib substitute for MATLAB's PATCH
    ('FaceColor','interp')/per-vertex FaceVertexCData, since
    Poly3DCollection only supports one flat colour per face.
    """
    return vertex_values[faces].mean(axis=1)
