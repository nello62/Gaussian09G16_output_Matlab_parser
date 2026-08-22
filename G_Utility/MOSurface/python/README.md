# MOSurface (Python port)

Python port of the MATLAB `G_Utility/MOSurface/` collection: real-space
molecular-orbital, electron-density, and electrostatic-potential
isosurfaces/contour maps evaluated directly from `.fchk` basis-set data,
plus load-and-replot from a saved Gaussian `.cube` file.

This is a standalone companion, mirroring the MATLAB original's own
status: **not** part of the `G09_*`/`G16_*` core toolbox, and **not**
folded into the tracked, published `G16parser` package. It lives here,
untracked, until a decision is made on how (or whether) to fold MO
visualization into the toolbox/paper.

## Naming philosophy

Same convention as `G16parser`: public functions are prefixed `g16_`
(there is only one Python package, so no `g09_`/`g16_` version split is
needed here either — see the manual's note on `g16_read_input` for why),
internal helpers are prefixed with a single leading underscore and are
not part of the public API.

## Public functions

| Function | MATLAB original |
|---|---|
| `g16_draw_mo_surface(data, mo_index, ...)` | `G_draw_mo_surface.m` |
| `g16_draw_density_surface(data, ...)` | `G_draw_density_surface.m` |
| `g16_draw_esp_surface(data, ...)` | `G_draw_esp_surface.m` |
| `g16_draw_cube_surface(cubefile, ...)` | `G_draw_cube_surface.m` |

All four accept the same options as their MATLAB counterparts, in
`snake_case` (e.g. `'IsoValue'` → `isovalue`, `'ShowMolecule'` →
`show_molecule`). `data` is the `Struct` returned by `G16parser`'s
`g16_fchk_read` — this is the one real dependency on `G16parser` (for
the input struct's shape and for the `g16_draw_molecule` CPK overlay),
mirroring the MATLAB original's own dependency on `G09_fchk_read`/
`G16_fchk_read` and `G09_draw_molecule`/`G16_draw_molecule`.

## Installation

```bash
pip install -r requirements.txt
# G16parser itself, if not already installed:
pip install -e ../../../G16parser
```

## Validation performed

The numerically critical parts were cross-checked directly against the
MATLAB original (not just internally self-consistent):

- `eval_mo_on_grid`/`eval_density_on_grid` (all S/P/SP/D/F/G shells,
  pure and Cartesian): agree with `g_eval_mo_on_grid.m`/
  `g_eval_density_on_grid.m` to ~1e-11 relative difference on a real
  `.fchk` (4-NTP, 218 basis functions), for HOMO/LUMO/HOMO-2/the
  highest-numbered MO and the total density.
- Pure D/F/G self-normalization: synthetic single-shell self-overlap
  integrals (fine-grid Riemann sum) all converge to 1.0000 for every
  component of every supported shell type.
- `eval_esp_analytic` (McMurchie-Davidson): agrees with
  `g_eval_esp_analytic.m` to ~4e-13 relative difference on the same
  real `.fchk`.
- Cube file read/write: round-trips exactly (to `%.5E` text precision)
  against a brute-force nested-loop reference, and cross-loads real
  cube files this port itself wrote.
- All four `g16_draw_*_surface` functions were run end-to-end on real
  `.fchk` files (4-NTP: closed-shell, S/P/SP/D shells; CH3 radical:
  UKS spin density; a field-on/field-off pair: density/ESP
  difference), producing chemically sensible pictures (e.g. the CH3
  radical spin-density plot shows the expected SOMO lobe perpendicular
  to the molecular plane; ESP-on-density coloring puts the
  electron-rich region at the expected heteroatom).

## Known differences from the MATLAB original (intentional, documented)

- **Isosurface extraction**: `skimage.measure.marching_cubes` in place
  of MATLAB's `isosurface`/`isonormals` (marching_cubes already returns
  per-vertex normals from the volume gradient, so no separate
  isonormals-equivalent step is needed).
- **Mesh decimation** (`max_vertices`, used by `g16_draw_esp_surface`
  and `g16_draw_cube_surface`): a simple vertex-clustering
  simplification (bin vertices on a coarse voxel grid, merge to
  centroids, remap faces) instead of MATLAB's `reducepatch`. A
  different, simpler algorithm, not a numerical-fidelity concern (ESP/
  ColorBy fields are smooth, so any reasonably-uniform vertex subset
  looks fine).
- **Shading**: matplotlib's `Poly3DCollection` has no true per-vertex
  Gouraud shading like MATLAB's `PATCH`; `shade=True` (simple per-face
  directional shading) is used for flat-coloured lobes, and a
  per-**face** colour (averaged from its 3 vertices) is used in place
  of MATLAB's `FaceColor,'interp'` per-vertex colouring for ESP/
  `color_by` surfaces.
- **`g16_draw_cube_surface`'s `color_by` interpolation** uses
  `scipy.interpolate.RegularGridInterpolator` with linear extrapolation
  at the domain edge (`fill_value=None`), rather than MATLAB
  `interp3`'s NaN-outside-domain behaviour: a marching-cubes vertex can
  land a floating-point hair outside the second cube's exact bounds,
  and since such a vertex is by construction already at the boundary,
  extrapolating there is effectively identical to the true boundary
  value.

## Performance characteristics (read before using on a large system)

Pure-Python/NumPy evaluation is markedly slower than MATLAB's compiled
engine for this workload, in two specific ways worth knowing about
before reaching for the defaults on a large molecule:

1. **Auto-expanding `padding`** (the default, `padding_explicit=False`):
   each retry re-evaluates the *entire* grid at a larger padding (up to
   4 attempts, ×1.6 each). In MATLAB this is cheap enough not to matter;
   in Python, on a real 218-basis-function molecule, one retry
   sequence alone was measured to take several minutes where an
   explicit, sufficient `padding` (no retries) took ~4 seconds.
   **Recommendation: pass an explicit `padding` (and
   `padding_explicit=True`) once you know a value that isn't clipped,
   rather than relying on the auto-expand default for routine use.**
2. **`g16_draw_esp_surface`'s analytic ESP evaluation** has a large
   *fixed* per-call cost dominated by the number of shell **pairs**
   (measured: ~11s fixed overhead for a 104-subshell basis, largely
   independent of how many points are evaluated, plus a smaller
   per-point marginal cost). This fixed cost is perfectly reasonable
   for the normal use case (a few hundred to a few thousand *decimated
   surface vertices*, via `max_vertices`), but makes `save_cube` +
   analytic ESP over a **full volumetric grid** (potentially hundreds
   of thousands of points) impractical at realistic resolution --
   unlike MATLAB, where the same feature is merely slow ("can be
   genuinely slow for a large basis", per `G_draw_esp_surface.m`'s own
   docstring) rather than impractical. **Recommendation: for
   `save_cube` on `g16_draw_esp_surface`, use a deliberately coarse
   `cube_spacing`/small `cube_padding`, or `esp_method='numeric'`
   (faster, less accurate -- see the MATLAB docstring's own "Accuracy"
   notes, which apply unchanged here) if a full-resolution ESP cube is
   genuinely needed.**

Optimizing the analytic ESP evaluator's inner shell-pair loop (e.g.
vectorizing across shell pairs instead of a nested Python loop) would
address point 2 but was out of scope for this initial port, which
prioritized exact numerical fidelity (validated to ~1e-13 relative
difference against MATLAB) over performance.
