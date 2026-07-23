# Gaussian09G16_output_Matlab_parser

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21302268.svg)](https://doi.org/10.5281/zenodo.21302268)
[![G16parser tests](https://github.com/nello62/Gaussian09G16_output_Matlab_parser/actions/workflows/python-tests.yml/badge.svg)](https://github.com/nello62/Gaussian09G16_output_Matlab_parser/actions/workflows/python-tests.yml)

MATLAB toolbox for parsing and visualising **Gaussian 09** and **Gaussian 16**
output files (`.out` / `.log` / `.fchk`) — molecular structure, energies,
charges, dipole moment and polarisability, vibrational normal modes, IR/Raman
spectra, orbital energies, and more.

The toolbox ships as two parallel, independent packages — **`G09/`** and
**`G16/`** — because the two Gaussian versions differ in output formatting.
Function names, signatures, and output struct fields are kept identical
between the two (`G09_xxx` / `G16_xxx`), so switching between a G09 and a G16
project only means changing which folder is on your MATLAB path.

Because the two formats differ, every data-extraction function checks the
file it reads against the Gaussian version it expects, and prints a
non-blocking warning if they do not match (e.g. calling `G16_energy` on a
Gaussian 09 file) — using the wrong toolbox can otherwise silently misparse
the file (e.g. polarisability coming back as `NaN`) rather than raising an
obvious error.

## Installation

```matlab
addpath('/path/to/Gaussian09G16_output_Matlab_parser/G09')   % for Gaussian 09 files
addpath('/path/to/Gaussian09G16_output_Matlab_parser/G16')   % for Gaussian 16 files
```

Add only the folder matching your files' Gaussian version (or both — the
`G09_*`/`G16_*` name prefixes never collide).

## Quick start

```matlab
mol = G16_structure('molecule.out');           % geometry
G16_draw_molecule(mol, 'ShowAxes', true);       % 3D CPK render

ch = G16_charges('molecule.out', 'ShowDipole', true);   % Mulliken/APT charges + dipole arrow

oe = G16_orbital_energies('molecule.out');      % HOMO/LUMO and full MO spectrum
G16_draw_orbital(oe);                           % energy-level diagram with HOMO-LUMO gap

nm = G16_nmodes('molecule.out');
G16_modeViewer('molecule.out');                 % interactive vibrational-mode browser

T = G16_read_all('molecule.out');               % everything above in one call
```

Every extraction function also prints a formatted summary to the command
window, and returns a struct so results can be used programmatically.

`read_all` reads the output file from disk only **once** and reuses the
parsed lines across every sub-function (via their shared `'Lines'`
parameter), instead of each one re-opening and re-parsing the file.

## Function reference

Run `G09_list()` / `G16_list()` at any time for the exact, up-to-date list
installed on your machine (returns a `table`, filterable by description).

### Data extraction

| Function | Description |
|---|---|
| `G_read_input` *(shared, no G09_/G16_ prefix)* | Reads a Gaussian **input** file (`.gjf`/`.com`/`.in`): link0 commands, route, title, charge/multiplicity, and starting geometry — output struct is compatible with `draw_molecule`/`get_bond_length` |
| `structure` | Molecular geometry (symbols, xyz, atom count) |
| `energy` | SCF energy and thermochemistry (ZPE, H, G, S, ...) |
| `charges` | Mulliken/APT atomic charges, optional dipole moment overlay |
| `dipole_polar` | Dipole moment and (static/dynamic) polarisability |
| `nmodes` | Vibrational normal-mode displacement vectors |
| `spectra` | IR and Raman spectra (stick + broadened) |
| `orbital_energies` | HOMO/LUMO and the full occupied/virtual MO spectrum |
| `convergence` | Geometry-optimisation convergence criteria per step |
| `charge_mult` | Molecular charge and spin multiplicity |
| `route` | Route section string |
| `get_bond_length` | Bond-length table from covalent radii |
| `gaussian_version` | Detects the Gaussian version/revision that produced the file (works on `.fchk` via a sibling `.log`/`.out`) |
| `read_all` | Runs the full extraction set in one call, reading the file only once |
| `batch_read_all` | Runs `read_all` (+ `orbital_energies`) over every `.log`/`.out` file in a folder and aggregates the key results into one summary table; per-file failures are recorded, not fatal |
| `fchk_read` *(G09 only)* | Reads a Gaussian formatted checkpoint (`.fchk`) file |
| `charges_fchk` *(G09 only)* | Visualises charges from a `fchk_read` struct |
| `restart` | Generates a `.gjf` restart input file from an existing output file |
| `hyperpolar` *(G16 only)* | Dipole hyperpolarisability (Beta) |
| `tddft` *(G16 only)* | TD-DFT excited states |

### Visualisation

| Function | Description |
|---|---|
| `draw_molecule` | 3D CPK ball-and-stick render, with optional Cartesian axes indicator; auto-detects double/triple bonds for C-C, C-N, C-O |
| `draw_mode` | 3D structure with a vibrational mode's displacement arrows |
| `draw_orbital` | Orbital energy-level diagram, HOMO-LUMO transition arrow + gap |
| `animate_mode` | Exports an MP4 animation of a vibrational mode (oscillating structure, GaussView-style) |
| `modeViewer` | Interactive mode selector/browser window, sortable by mode number, IR, or Raman intensity, with an "Animate mode (MP4)..." button |

### Utility

| Function | Description |
|---|---|
| `list` | Lists every function in the toolbox with its one-line description |
| `write_report` | Writes a formatted text report (.txt) from a `read_all` struct |
| `read_lines` *(G09 only)* | Shared file-reading helper used internally |
| `check_gaussian_match` | Internal: warns if a file looks like the other Gaussian version (called automatically by every reading function) |

All Name-Value options, output struct fields, and examples are documented in
each function's own help text (`help G16_charges`, etc.).

## Python port

A Python 3 port of the G16 toolbox lives in [`G16parser/`](G16parser/) —
data extraction, static matplotlib plots, and an interactive Tkinter
vibrational-mode viewer (`g16_mode_viewer`). Installable as a standard
package:

```bash
cd G16parser
pip install -e .
```

```python
import G16parser as g16

mol = g16.g16_structure('molecule.out')
g16.g16_draw_molecule(mol, show_axes=True)

T = g16.g16_read_all('molecule.out')   # everything in one call, single file read
```

See [`G16parser/README.md`](G16parser/README.md) for the full function
reference and [`G16parser/example.py`](G16parser/example.py) for a runnable
end-to-end example.

## License

See [`LICENSE`](LICENSE).
