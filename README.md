# Gaussian09G16_output_Matlab_parser

MATLAB toolbox for parsing and visualising **Gaussian 09** and **Gaussian 16**
output files (`.out` / `.log` / `.fchk`) — molecular structure, energies,
charges, dipole moment and polarisability, vibrational normal modes, IR/Raman
spectra, orbital energies, and more.

The toolbox ships as two parallel, independent packages — **`G09/`** and
**`G16/`** — because the two Gaussian versions differ in output formatting.
Function names, signatures, and output struct fields are kept identical
between the two (`G09_xxx` / `G16_xxx`), so switching between a G09 and a G16
project only means changing which folder is on your MATLAB path.

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
```

Every extraction function also prints a formatted summary to the command
window, and returns a struct so results can be used programmatically.

## Function reference

Run `G09_list()` / `G16_list()` at any time for the exact, up-to-date list
installed on your machine (returns a `table`, filterable by description).

### Data extraction

| Function | Description |
|---|---|
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
| `fchk_read` *(G09 only)* | Reads a Gaussian formatted checkpoint (`.fchk`) file |
| `charges_fchk` *(G09 only)* | Visualises charges from a `fchk_read` struct |
| `restart` *(G09 only)* | Generates a `.gjf` restart input file |
| `hyperpolar` *(G16 only)* | Dipole hyperpolarisability (Beta) |
| `tddft` *(G16 only)* | TD-DFT excited states |

### Visualisation

| Function | Description |
|---|---|
| `draw_molecule` | 3D CPK ball-and-stick render, with optional Cartesian axes indicator |
| `draw_mode` | 3D structure with a vibrational mode's displacement arrows |
| `draw_orbital` | Orbital energy-level diagram, HOMO-LUMO transition arrow + gap |
| `modeViewer` | Interactive mode selector/browser window, sortable by mode number, IR, or Raman intensity |

### Utility

| Function | Description |
|---|---|
| `list` | Lists every function in the toolbox with its one-line description |
| `read_lines` *(G09 only)* | Shared file-reading helper used internally |

All Name-Value options, output struct fields, and examples are documented in
each function's own help text (`help G16_charges`, etc.).

## License

See [`LICENSE`](LICENSE).
