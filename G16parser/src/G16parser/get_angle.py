import numpy as np
import pandas as pd


def _coords(mol, fname):
    """(xyz [N x 3] ndarray, symbols list) from a Struct with .xyz (and
    optionally .symbols) or from an N x 3 array."""
    if hasattr(mol, "xyz"):
        xyz = np.asarray(mol.xyz, dtype=float)
        sym = list(getattr(mol, "symbols", []))
        if len(sym) != xyz.shape[0]:
            sym = [""] * xyz.shape[0]
    else:
        xyz = np.asarray(mol, dtype=float)
        sym = [""] * xyz.shape[0]
    if xyz.ndim != 2 or xyz.shape[1] != 3:
        raise ValueError(f"{fname}: coordinates must be an N x 3 array.")
    return xyz, sym


def _indices(atoms, ncol, natoms, fname):
    a = np.atleast_2d(np.asarray(atoms))
    if a.shape[1] != ncol:
        raise ValueError(f"{fname}: atoms must have {ncol} columns.")
    if np.any(a != np.round(a)) or a.min() < 1 or a.max() > natoms:
        raise ValueError(f"{fname}: atom indices must be integers between 1 and {natoms}.")
    return a.astype(int)


def g16_get_angle(mol, atoms, save_as="", verbose=False):
    """Bond angles i-j-k (degrees) from a molecular structure.

    Parameters
    ----------
    mol : Struct with .xyz (and optionally .symbols), as returned by
        g16_structure / g16_read_input, or an N x 3 coordinate array
    atoms : [M x 3] (or a single triple) of 1-based atom indices; each row
        defines the angle atoms[m,0]-atoms[m,1]-atoms[m,2] (vertex in the
        middle)
    save_as : str, optional — .xlsx or .csv file for the table
    verbose : bool — print the table (default False)

    Returns
    -------
    ang_table : pandas.DataFrame with columns Atom1, Sym1, Atom2, Sym2,
        Atom3, Sym3, Angle_deg (1-based indices, as in the MATLAB output)
    ang : numpy.ndarray [M] of angles in degrees (0..180), computed as
        atan2(|u x v|, u . v) with u = r1 - r2, v = r3 - r2
    """
    xyz, sym = _coords(mol, "g16_get_angle")
    a = _indices(atoms, 3, xyz.shape[0], "g16_get_angle")
    ang = np.empty(len(a))
    for m, (i, j, k) in enumerate(a):
        u = xyz[i - 1] - xyz[j - 1]
        v = xyz[k - 1] - xyz[j - 1]
        ang[m] = np.degrees(np.arctan2(np.linalg.norm(np.cross(u, v)), np.dot(u, v)))
    rows = [(i, sym[i - 1], j, sym[j - 1], k, sym[k - 1], x) for (i, j, k), x in zip(a, ang)]
    ang_table = pd.DataFrame(rows, columns=["Atom1", "Sym1", "Atom2", "Sym2",
                                            "Atom3", "Sym3", "Angle_deg"])
    if verbose:
        print(f"\n{'angle':16s} {'deg':>10s}")
        for r in rows:
            print(f"{r[1]}{r[0]}-{r[3]}{r[2]}-{r[5]}{r[4]}".ljust(16) + f" {r[6]:10.3f}")
        print()
    if save_as:
        if save_as.lower().endswith(".csv"):
            ang_table.to_csv(save_as, index=False)
        else:
            ang_table.to_excel(save_as, index=False)
        print(f"Table saved to: {save_as}")
    return ang_table, ang
