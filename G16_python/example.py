"""Basic usage example for the G16_python toolbox.

Run with:
    python3 example.py path/to/molecule.out
"""
import os
import sys

# Make "import G16_python" work when this file is run directly
# (python3 example.py ...), not just as a module (-m G16_python.example).
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib.pyplot as plt

import G16_python as g16


def main(filename):
    # Geometry
    mol = g16.g16_structure(filename)
    print(f"{mol.Natoms} atoms, formula: {'-'.join(sorted(set(mol.symbols)))}")

    # SCF energy and thermochemistry
    en = g16.g16_energy(filename)
    print(f"SCF energy: {en.SCF:.6f} Ha")

    # Atomic charges + dipole moment (plot=False: data only, no figure)
    ch = g16.g16_charges(filename, plot=False, show_dipole=True)
    print(f"Dipole moment: {ch.dipole_Debye:.3f} D")

    # HOMO/LUMO
    oe = g16.g16_orbital_energies(filename)
    print(f"HOMO-LUMO gap: {oe.gap_eV:.3f} eV")

    # --- Figures -----------------------------------------------------
    g16.g16_draw_molecule(mol, show_axes=True)
    g16.g16_draw_orbital(oe)

    # Tip: g16.g16_read_all(filename) runs everything above (and more)
    # in a single call, reading the file from disk only once.

    plt.show()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "molecule.out")
