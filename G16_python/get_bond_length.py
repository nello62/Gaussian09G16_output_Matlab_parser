import warnings

import numpy as np
import pandas as pd

# Note: this table is intentionally separate from _common.COVALENT_RADII
# (used by draw_molecule) — same source data, but a different unknown-
# element fallback (1.50 A + warning here, vs. a silent 0.80 A default in
# draw_molecule), mirroring the two independent tables in the MATLAB code.
_ELEMENTS = ["H", "C", "N", "O", "F", "S", "P", "Cl", "Br", "I",
             "Si", "B", "Na", "K", "Mg", "Ca",
             "Au", "Ag", "Cu", "Zn", "Cd", "Pt", "Pd", "Ni", "Fe", "Co"]
_RADII = [0.31, 0.76, 0.71, 0.66, 0.57, 1.05, 1.07, 1.02, 1.20, 1.39,
          1.11, 0.84, 1.66, 2.03, 1.41, 1.76,
          1.36, 1.45, 1.32, 1.22, 1.44, 1.36, 1.39, 1.24, 1.32, 1.26]
_COVALENT_RADII = dict(zip(_ELEMENTS, _RADII))


def _lookup_radius(sym):
    if sym in _COVALENT_RADII:
        return _COVALENT_RADII[sym]
    warnings.warn(
        f'Elemento "{sym}" non presente nella tabella dei raggi covalenti: '
        f"uso raggio di default 1.50 A."
    )
    return 1.50


def g16_get_bond_length(mol, tolerance=1.15, include_h=True, sort_by="distance", save_as=""):
    """Builds the bond-length table of a molecule using a covalent-radii
    geometric criterion (not derived from Gaussian's own connectivity).

    Parameters
    ----------
    mol : Struct — as returned by g16_structure (needs .symbols, .xyz)
    tolerance : float — bonded if dist <= tolerance * (r1 + r2) (default 1.15)
    include_h : bool — include bonds involving hydrogens (default True)
    sort_by : 'distance' (default) | 'atom'
    save_as : str, optional — path (.xlsx or .csv) to save the table to

    Returns
    -------
    bond_table : pandas.DataFrame with columns Atom1, Sym1, Atom2, Sym2,
        Distance_Ang (Atom1/Atom2 are 1-based, matching the MATLAB output).
    """
    symbols = list(mol.symbols)
    xyz = np.asarray(mol.xyz)
    natoms = len(symbols)

    if xyz.shape[0] != natoms:
        raise ValueError(
            f"g16_get_bond_length: il numero di simboli ({natoms}) non corrisponde "
            f"al numero di righe di xyz ({xyz.shape[0]})."
        )

    rows = []
    for i in range(natoms - 1):
        si = symbols[i]
        if not include_h and si.upper() == "H":
            continue
        ri = _lookup_radius(si)
        for j in range(i + 1, natoms):
            sj = symbols[j]
            if not include_h and sj.upper() == "H":
                continue
            rj = _lookup_radius(sj)
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            if d <= tolerance * (ri + rj):
                rows.append((i + 1, si, j + 1, sj, d))

    bond_table = pd.DataFrame(rows, columns=["Atom1", "Sym1", "Atom2", "Sym2", "Distance_Ang"])

    if sort_by.lower() == "distance":
        bond_table = bond_table.sort_values("Distance_Ang").reset_index(drop=True)
    elif sort_by.lower() == "atom":
        bond_table = bond_table.sort_values(["Atom1", "Atom2"]).reset_index(drop=True)

    if save_as:
        if save_as.lower().endswith(".csv"):
            bond_table.to_csv(save_as, index=False)
        else:
            bond_table.to_excel(save_as, index=False)
        print(f"Tabella dei legami salvata in: {save_as}")

    return bond_table
