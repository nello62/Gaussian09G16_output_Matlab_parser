# G16parser

Python 3 port of the `G16/` MATLAB toolbox — parses and visualises Gaussian 16
`.out`/`.log`/`.fchk` files. Data-extraction functions, static matplotlib
plots, and an interactive Tkinter vibrational-mode viewer.

## Requirements

```
numpy
pandas
matplotlib
```

## Install / use

Installable as a standard Python package (editable install, so local edits
take effect immediately without reinstalling):

```
cd G16parser
pip install -e .
```

Then, from any folder/script:

```python
import G16parser as g16

mol = g16.g16_structure('molecule.out')
g16.g16_draw_molecule(mol, show_axes=True)

ginp = g16.g16_read_input('molecule.gjf')   # also .com / .in
g16.g16_draw_molecule(ginp, show_axes=True)  # reuses the same struct shape

ch = g16.g16_charges('molecule.out', show_dipole=True)

oe = g16.g16_orbital_energies('molecule.out')
g16.g16_draw_orbital(oe)

T = g16.g16_read_all('molecule.out')   # everything in one call, single file read
g16.g16_write_report(T)                # -> molecule_report.txt

g16.g16_restart('molecule.out')        # -> molecule_restarted.gjf

summary = g16.g16_batch_read_all('results/', write_reports=True)  # one row per file
```

Every function returns a `Struct` (attribute access, e.g. `mol.xyz`,
`ch.dipole_Debye`) — the Python equivalent of the MATLAB struct outputs.

See [`example.py`](example.py) for a short, runnable end-to-end example:

```
python3 example.py path/to/molecule.out
```

## Running tests

```bash
cd G16parser
pip install -e ".[test]"
pytest
```

Most tests need one real Gaussian 16 `.out`/`.log` file in
`tests/fixtures/` to run against (see
[`tests/fixtures/README.md`](tests/fixtures/README.md)); without one they
skip cleanly rather than failing.

## Function reference

| Function | Description |
|---|---|
| `g16_read_input` | Reads a Gaussian **input** file (`.gjf`/`.com`/`.in`): link0 commands, route, title, charge/multiplicity, and starting geometry — output `Struct` is compatible with `g16_draw_molecule`/`g16_get_bond_length` |
| `g16_structure` | Molecular geometry (symbols, xyz as np.ndarray, atom count) |
| `g16_energy` | SCF energy and thermochemistry |
| `g16_convergence` | Geometry-optimisation convergence criteria per step (`plot=True` for a matplotlib figure) |
| `g16_charges` | Mulliken/APT atomic charges, optional dipole overlay and 3D plot |
| `g16_dipole_polar` | Dipole moment and (static/dynamic) polarisability |
| `g16_nmodes` | Vibrational normal-mode displacement vectors |
| `g16_spectra` | IR/Raman spectra (stick + Lorentzian-broadened continuum) |
| `g16_orbital_energies` | HOMO/LUMO and the full occupied/virtual MO spectrum |
| `g16_charge_mult` | Molecular charge and spin multiplicity |
| `g16_route` | Route section string |
| `g16_get_bond_length` | Bond-length table (pandas DataFrame) from covalent radii |
| `g16_gaussian_version` | Detects the Gaussian version/revision (works on `.fchk` via a sibling `.log`/`.out`) |
| `g16_hyperpolar` | Dipole hyperpolarisability (Beta) |
| `g16_tddft` | TD-DFT excited states |
| `g16_read_all` | Runs the full extraction set in one call, reading the file only once |
| `g16_restart` | Generates a `.gjf` restart input file from an existing output file |
| `g16_batch_read_all` | Runs `g16_read_all` (+ `g16_orbital_energies`) over every `.log`/`.out` file in a folder and aggregates the key results into one summary DataFrame; per-file failures are recorded, not fatal |
| `g16_write_report` | Writes a formatted text report (.txt) from a `g16_read_all` Struct |
| `g16_draw_molecule` | 3D CPK ball-and-stick render (matplotlib 3D); auto-detects double/triple bonds for C-C, C-N, C-O |
| `g16_draw_mode` | 3D structure with a vibrational mode's displacement arrows |
| `g16_draw_orbital` | Orbital energy-level diagram with HOMO-LUMO gap arrow |
| `g16_animate_mode` | Exports an MP4 animation of a vibrational mode (requires `ffmpeg`) |
| `g16_mode_viewer` | Interactive Tkinter window to browse/render vibrational modes, sortable by mode number, IR, or Raman intensity, with an "Animate mode (MP4)..." button |
| `g16_list` | Lists every function in the toolbox with its one-line description (pandas DataFrame) |

## Notes on the port

- Single-file-read optimisation: every extraction function accepts an
  optional `lines=` parameter (pre-read file lines) to avoid re-reading
  the file, mirroring the MATLAB toolbox's `'Lines'` parameter; `g16_read_all`
  uses this to read the file exactly once.
- 3D plots use matplotlib's `mplot3d`; there is no MATLAB-style
  `rotate3d` — use the interactive backend's normal mouse/toolbar controls.
- The package forces the `TkAgg` matplotlib backend (if Tk is available and
  pyplot hasn't been imported yet) instead of the native macOS backend: on
  some older matplotlib releases, interactively rotating a 3D plot with the
  native Cocoa backend can segfault the whole Python process. If you need a
  different backend, call `matplotlib.use(...)` yourself *before* importing
  `G16parser`.
- A bug discovered in the original `G16_tddft.m` (fixed 2026-07-13): MATLAB's
  `regexp` drops a trailing optional `<S**2>=` capture group entirely (a
  MATLAB-specific regex-engine quirk — confirmed with `'tokens','once'`,
  plain `'tokens'`, and `'names'`), so `td.S2` used to always be `NaN` in
  MATLAB even when the tag was present. Fixed in MATLAB by extracting
  `<S**2>=` with its own independent `regexp` call instead of a trailing
  optional group. This Python port's `re`-based parser never had the bug
  (Python's `re` handles trailing optional groups correctly) and has always
  extracted `S2` correctly.
- `g16_mode_viewer` is a Tkinter rewrite of `G16_modeViewer.m`'s `uifigure`
  GUI (same controls: mode selector, order-by, per-option redraw, title,
  save-as PDF/EPS/JPEG). One simplification: "target figure" always means
  the most recently drawn mode window, not whichever window the user last
  clicked on (MATLAB's `groot().CurrentFigure` has no simple Tkinter
  equivalent). In an IPython/Spyder console with the Tk graphics backend
  active, the viewer detects the running event loop and returns to the
  prompt immediately instead of blocking; run as a plain script, it blocks
  until the viewer window is closed.
- **Known issue:** `g16_mode_viewer` can segfault when run from *inside*
  Spyder's IPython console (crash inside ipykernel's own periodic Tk
  event-loop pump, not in this toolbox's code) — confirmed working
  correctly when run from a plain terminal/script instead:
  `python3 -c "import G16parser as g16; g16.g16_mode_viewer('file.out')"`.
  All other functions (including static plots) work fine inside Spyder.
- `g16_animate_mode` requires **`ffmpeg`** installed and on `PATH` —
  matplotlib does not bundle a video encoder itself:
  `brew install ffmpeg` (macOS) or `sudo apt install ffmpeg` (Ubuntu/Debian).
  Without it, the call fails with `FileNotFoundError: ... 'ffmpeg'` at the
  final encoding step (frame generation itself does not require ffmpeg).
  By default the animation uses matplotlib's default 3D view; pass
  `view=(ax.azim, ax.elev)` to start from a figure's current orientation
  instead (`g16_mode_viewer`'s "Animate mode (MP4)..." button does this
  automatically, using whichever mode figure is currently displayed).
  Bonds are computed once from the equilibrium geometry (via
  `g16_get_bond_length`) and kept fixed across all frames, instead of
  being re-detected from instantaneous distances every frame — otherwise
  bonds flicker in and out as atoms oscillate past the `bond_tol`
  threshold.
- **Bond order (single/double/triple):** `g16_draw_molecule` estimates
  bond order purely from bond length for C-C, C-N, and C-O pairs (any
  other element pair is always drawn as a single bond), rendered as
  1/2/3 parallel lines — a geometric estimate, not Gaussian's own
  bond-order analysis (e.g. Wiberg/NBO indices). The `bond_list`
  parameter accepts an optional 3rd column with a pre-computed order,
  used by `g16_animate_mode` to keep both bond topology and order fixed
  across all frames. The C-C double/single threshold is set to 1.36 Å
  rather than the generic reference-length midpoint, so symmetric
  aromatic rings (C-C ≈1.39-1.40 Å, no real length alternation) render
  as all-single instead of all-double.
