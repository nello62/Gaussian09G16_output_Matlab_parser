"""Isosurface extraction/decimation and CPK molecule rendering shared by
all four g16_draw_*_surface functions.

Isosurface extraction uses skimage.measure.marching_cubes as the
matplotlib-ecosystem equivalent of MATLAB's built-in isosurface/
isonormals (marching_cubes already returns per-vertex normals from the
volume gradient, so no separate isonormals-equivalent step is needed).

Mesh decimation (MATLAB's REDUCEPATCH, used by G_draw_esp_surface's
'MaxVertices' and G_draw_cube_surface's own) has no scikit-image
equivalent; a simple, clearly-documented vertex-clustering simplification
is used instead (bin vertices on a coarse voxel grid, merge each bin to
its centroid, remap faces, drop degenerate ones) -- this is a genuine,
intentional implementation difference from MATLAB (a different,
simpler algorithm, not a numerical-fidelity concern: ESP/ColorBy fields
are smooth, so any reasonably-uniform vertex subset gives a visually
smooth result).
"""

import numpy as np
from skimage import measure


def marching_cubes_ang(X_ang, Y_ang, Z_ang, V, isovalue):
    """Extracts the isosurface of V at isovalue, in Angstrom coordinates
    (X_ang/Y_ang/Z_ang are 1D axis vectors, V has meshgrid convention
    V[iy,ix,iz] matching np.meshgrid(X_ang,Y_ang,Z_ang) default 'xy'
    indexing).

    Returns (vertices [Nv,3], faces [Nf,3], normals [Nv,3]) in Angstrom,
    or (None, None, None) if no isosurface is found at this level (V does
    not span isovalue).
    """
    vmin, vmax = V.min(), V.max()
    if not (vmin < isovalue < vmax):
        return None, None, None

    dx = X_ang[1] - X_ang[0] if len(X_ang) > 1 else 1.0
    dy = Y_ang[1] - Y_ang[0] if len(Y_ang) > 1 else 1.0
    dz = Z_ang[1] - Z_ang[0] if len(Z_ang) > 1 else 1.0

    # marching_cubes expects volume[i,j,k] with spacing along (i,j,k);
    # V is (Ny,Nx,Nz) -- pass directly, spacing=(dy,dx,dz), then swap the
    # first two returned vertex columns to get (x,y,z) order.
    verts, faces, normals, _values = measure.marching_cubes(V, level=isovalue, spacing=(dy, dx, dz))
    verts_xyz = np.column_stack([verts[:, 1] + X_ang[0], verts[:, 0] + Y_ang[0], verts[:, 2] + Z_ang[0]])
    normals_xyz = np.column_stack([normals[:, 1], normals[:, 0], normals[:, 2]])
    return verts_xyz, faces, normals_xyz


def decimate_mesh(vertices, faces, max_vertices, normals=None, values=None):
    """Simple vertex-clustering mesh decimation to (approximately) at
    most max_vertices vertices -- see this module's docstring for why
    this differs from MATLAB's REDUCEPATCH. Optionally carries per-vertex
    normals and/or a scalar value array through the same merge (averaged
    per cluster, normals renormalized).

    Returns a dict with keys 'vertices', 'faces', and (if given)
    'normals'/'values'.
    """
    Nv = vertices.shape[0]
    if Nv <= max_vertices:
        out = {"vertices": vertices, "faces": faces}
        if normals is not None:
            out["normals"] = normals
        if values is not None:
            out["values"] = values
        return out

    lo = vertices.min(axis=0)
    hi = vertices.max(axis=0)
    extent = np.maximum(hi - lo, 1e-9)
    # Aim for ~max_vertices occupied bins by choosing a per-axis bin
    # count proportional to the bounding box's aspect ratio.
    n_bins_total = max(max_vertices, 8)
    bins_per_axis = max(2, int(round(n_bins_total ** (1 / 3))))

    idx = np.floor((vertices - lo) / extent * (bins_per_axis - 1e-9)).astype(int)
    idx = np.clip(idx, 0, bins_per_axis - 1)
    keys = (idx[:, 0] * bins_per_axis + idx[:, 1]) * bins_per_axis + idx[:, 2]

    uniq_keys, inverse = np.unique(keys, return_inverse=True)
    n_new = len(uniq_keys)

    new_vertices = np.zeros((n_new, 3))
    counts = np.bincount(inverse, minlength=n_new)
    for d in range(3):
        new_vertices[:, d] = np.bincount(inverse, weights=vertices[:, d], minlength=n_new) / counts

    new_faces_raw = inverse[faces]
    degenerate = (new_faces_raw[:, 0] == new_faces_raw[:, 1]) | \
                 (new_faces_raw[:, 1] == new_faces_raw[:, 2]) | \
                 (new_faces_raw[:, 0] == new_faces_raw[:, 2])
    new_faces = new_faces_raw[~degenerate]

    out = {"vertices": new_vertices, "faces": new_faces}

    if normals is not None:
        new_normals = np.zeros((n_new, 3))
        for d in range(3):
            new_normals[:, d] = np.bincount(inverse, weights=normals[:, d], minlength=n_new) / counts
        norm_len = np.linalg.norm(new_normals, axis=1, keepdims=True)
        norm_len[norm_len == 0] = 1.0
        out["normals"] = new_normals / norm_len

    if values is not None:
        out["values"] = np.bincount(inverse, weights=values, minlength=n_new) / counts

    return out
