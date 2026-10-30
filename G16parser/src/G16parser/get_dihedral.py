import warnings

import numpy as np
import pandas as pd

from .get_angle import _coords, _indices


def g16_get_dihedral(mol, atoms, signed=True, save_as="", verbose=False):
    """Dihedral (torsion) angles i-j-k-l (degrees) from a molecular structure.

    Parameters
    ----------
    mol : Struct with .xyz (and optionally .symbols), as returned by
        g16_structure / g16_read_input, or an N x 3 coordinate array
    atoms : [M x 4] (or a single quadruple) of 1-based atom indices; each
        row defines the angle between the planes (1,2,3) and (2,3,4)
    signed : bool — True (default): IUPAC sign, range (-180, 180];
        False: absolute value, range [0, 180]
    save_as : str, optional — .xlsx or .csv file for the table
    verbose : bool — print the table (default False)

    Returns
    -------
    dih_table : pandas.DataFrame with columns Atom1, Sym1, ..., Atom4, Sym4,
        Dihedral_deg (1-based indices, as in the MATLAB output)
    phi : numpy.ndarray [M] of dihedral angles in degrees, from
        atan2((n1 x n2) . b2/|b2|, n1 . n2) with b1 = r2 - r1,
        b2 = r3 - r2, b3 = r4 - r3, n1 = b1 x b2, n2 = b2 x b3;
        NaN (with a warning) if three consecutive atoms are collinear
    """
    xyz, sym = _coords(mol, "g16_get_dihedral")
    a = _indices(atoms, 4, xyz.shape[0], "g16_get_dihedral")
    phi = np.empty(len(a))
    for m, idx in enumerate(a):
        r = xyz[idx - 1]
        b1, b2, b3 = r[1] - r[0], r[2] - r[1], r[3] - r[2]
        n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
        tol = 1e-8
        if (np.linalg.norm(n1) < tol * np.linalg.norm(b1) * np.linalg.norm(b2)
                or np.linalg.norm(n2) < tol * np.linalg.norm(b2) * np.linalg.norm(b3)):
            warnings.warn(f"Dihedral {'-'.join(map(str, idx))} undefined: "
                          "three consecutive atoms are collinear.")
            phi[m] = np.nan
            continue
        phi[m] = np.degrees(np.arctan2(np.dot(np.cross(n1, n2), b2 / np.linalg.norm(b2)),
                                       np.dot(n1, n2)))
    if not signed:
        phi = np.abs(phi)
    rows = [(i, sym[i - 1], j, sym[j - 1], k, sym[k - 1], l, sym[l - 1], x)
            for (i, j, k, l), x in zip(a, phi)]
    dih_table = pd.DataFrame(rows, columns=["Atom1", "Sym1", "Atom2", "Sym2", "Atom3",
                                            "Sym3", "Atom4", "Sym4", "Dihedral_deg"])
    if verbose:
        print(f"\n{'dihedral':22s} {'deg':>10s}")
        for r in rows:
            lab = f"{r[1]}{r[0]}-{r[3]}{r[2]}-{r[5]}{r[4]}-{r[7]}{r[6]}"
            print(lab.ljust(22) + f" {r[8]:10.3f}")
        print()
    if save_as:
        if save_as.lower().endswith(".csv"):
            dih_table.to_csv(save_as, index=False)
        else:
            dih_table.to_excel(save_as, index=False)
        print(f"Table saved to: {save_as}")
    return dih_table, phi
