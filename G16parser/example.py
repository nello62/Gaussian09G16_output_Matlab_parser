"""Basic usage example for the G16parser toolbox.

Requires G16parser to be installed (pip install -e . from this folder).

Run with:
    python3 example.py path/to/molecule.out
"""
import sys

# Import G16parser BEFORE matplotlib.pyplot: the package switches the
# matplotlib backend to TkAgg on import (see _common.py) to avoid a
# native-macOS-backend segfault during interactive 3D rotation, but it
# can only do that if pyplot has not been imported yet.
import G16parser as g16
import matplotlib.pyplot as plt


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
    # Mode
    nm=g16.g16_nmodes(filename) 
    
    # --- Figures -----------------------------------------------------
    g16.g16_draw_molecule(mol, show_axes=True)
    g16.g16_draw_orbital(oe)
    g16.g16_draw_mode(mol, nm, 60)
    
    # Tip: g16.g16_read_all(filename) runs everything above (and more)
    # in a single call, reading the file from disk only once.

    plt.show()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "molecule.out")
