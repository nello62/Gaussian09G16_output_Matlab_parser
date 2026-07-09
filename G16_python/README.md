# G16_python

Python 3 port of the `G16/` MATLAB toolbox — parses and visualises Gaussian 16
`.out`/`.log`/`.fchk` files. Data-extraction functions plus static
matplotlib plots; the interactive `G16_modeViewer` GUI is **not** ported.

## Requirements

```
numpy
pandas
matplotlib
```

## Install / use

No packaging yet — just make sure the parent folder is on `sys.path` (or
run from `G16_python`'s parent directory) and import directly:

```python
import G16_python as g16

mol = g16.g16_structure('molecule.out')
g16.g16_draw_molecule(mol, show_axes=True)

ch = g16.g16_charges('molecule.out', show_dipole=True)

oe = g16.g16_orbital_energies('molecule.out')
g16.g16_draw_orbital(oe)

T = g16.g16_read_all('molecule.out')   # everything in one call, single file read
```

Every function returns a `Struct` (attribute access, e.g. `mol.xyz`,
`ch.dipole_Debye`) — the Python equivalent of the MATLAB struct outputs.

See [`example.py`](example.py) for a short, runnable end-to-end example:

```
python3 example.py path/to/molecule.out
```

## Function reference

| Function | Description |
|---|---|
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
| `g16_draw_molecule` | 3D CPK ball-and-stick render (matplotlib 3D) |
| `g16_draw_mode` | 3D structure with a vibrational mode's displacement arrows |
| `g16_draw_orbital` | Orbital energy-level diagram with HOMO-LUMO gap arrow |

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
  `G16_python`.
- A discovered bug in the original `G16_tddft.m`: MATLAB's `regexp` drops
  the optional trailing `<S**2>=` capture group entirely (a MATLAB-specific
  regex-engine quirk), so `td.S2` is always `NaN` in MATLAB even when the
  tag is present in the file. This Python port's `re`-based parser does not
  have that bug and extracts `S2` correctly.
