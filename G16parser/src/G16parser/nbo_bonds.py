import re

import pandas as pd

from ._common import read_lines

_HDR_RE = re.compile(r"Natural Bond Orbitals \(Summary\)")
# Bonding (BD) NBO lines look like:
#     1. BD (   1) C   1 - H   2          1.99909    -0.50504
# Antibonding (BD*) lines are deliberately excluded: the pattern requires
# "BD" to be followed only by whitespace then "(", so "BD*(" cannot match.
_BD_RE = re.compile(
    r"^\s*\d+\.\s+BD\s*\(\s*\d+\)\s+([A-Za-z]+)\s*(\d+)\s*-\s*([A-Za-z]+)\s*(\d+)\s+([\d.]+)"
)


def g16_nbo_bonds(filename, section="last", lines=None):
    """Determines bond order (single/double/triple) from a Gaussian NBO
    (Natural Bond Orbital) analysis.

    Requires the source file to have been computed with the 'pop=nbo'
    Gaussian keyword (or similar, e.g. 'pop=(nbo,savenbo)'); without it,
    the output file simply does not contain NBO data and this function
    raises a ValueError -- bond order cannot be recovered after the fact
    from geometry/energy alone, the NBO analysis has to have actually run
    as part of the job.

    Unlike g16_get_bond_length (a purely geometric covalent-radius
    criterion) and g16_draw_molecule's bond-order heuristic (a bond
    length threshold for C-C/C-O pairs only, explicitly not derived from
    any real Gaussian bond-order analysis), this reads the actual NBO
    "Natural Bond Orbitals (Summary)" table and counts, for every atom
    pair, how many bonding (BD) natural bond orbitals NBO assigned to it:
    one BD orbital = single bond (sigma only), two = double bond (sigma +
    pi), three = triple bond (sigma + 2 pi). This reflects Gaussian's own
    NBO analysis, not a geometric guess.

    Note: this reads the "Natural Bond Orbitals (Summary)" table that
    'pop=nbo' always prints, not the separate "Wiberg bond index matrix"
    (which additionally requires 'pop=(nbo,bndidx)' and is not parsed by
    this function).

    Parameters
    ----------
    filename : str
    section : 'last' (default) | 'first' — which "Natural Bond Orbitals
        (Summary)" block to read, for files with more than one (e.g.
        multi-step opt+freq+polar jobs that repeat the NBO analysis at
        each step)
    lines : list[str], optional — pre-read file lines, to skip re-reading
        the file when it has already been read elsewhere.

    Returns
    -------
    bond_table : pandas.DataFrame with columns Atom1, Sym1, Atom2, Sym2
        (Atom1/Atom2 1-based, Atom1 < Atom2), BondOrder (int, count of BD
        orbitals for this pair), Occupancy (summed NBO occupancy, ~2 per
        BD orbital).

    Example
    -------
        bt = g16_nbo_bonds('CH4_NBO.LOG')
    """
    if lines is None:
        lines = read_lines(filename)

    hdr_idx = [i for i, ln in enumerate(lines) if _HDR_RE.search(ln)]
    if not hdr_idx:
        raise ValueError(
            f"g16_nbo_bonds: no NBO analysis found in {filename}. This requires "
            "the source Gaussian job to have been run with the 'pop=nbo' keyword "
            "(or similar) -- without it, bond-order data is not present in the file."
        )

    k0 = hdr_idx[-1] if section.lower() == "last" else hdr_idx[0]

    pairs = []
    for ln in lines[k0:]:
        if "Total Lewis" in ln:
            break
        m = _BD_RE.match(ln)
        if not m:
            continue
        raw_s1, raw_a1, raw_s2, raw_a2, occ_s = m.groups()
        raw_a1, raw_a2, occ = int(raw_a1), int(raw_a2), float(occ_s)
        if raw_a1 <= raw_a2:
            a1, s1, a2, s2 = raw_a1, raw_s1, raw_a2, raw_s2
        else:
            a1, s1, a2, s2 = raw_a2, raw_s2, raw_a1, raw_s1
        pairs.append((a1, s1, a2, s2, occ))

    if not pairs:
        raise ValueError(
            f"g16_nbo_bonds: NBO summary section found in {filename} but no "
            "bonding (BD) orbitals were read."
        )

    df = pd.DataFrame(pairs, columns=["Atom1", "Sym1", "Atom2", "Sym2", "_occ"])
    grouped = df.groupby(["Atom1", "Atom2"], sort=False)
    bond_table = grouped.agg(
        Sym1=("Sym1", "first"),
        Sym2=("Sym2", "first"),
        BondOrder=("_occ", "size"),
        Occupancy=("_occ", "sum"),
    ).reset_index()
    bond_table = bond_table[["Atom1", "Sym1", "Atom2", "Sym2", "BondOrder", "Occupancy"]]
    bond_table = bond_table.sort_values(["Atom1", "Atom2"]).reset_index(drop=True)

    print(f"\n-- g16_nbo_bonds: {filename} --")
    print(f"  {len(bond_table)} bonds found from NBO analysis (section: {section})")
    for _, row in bond_table.iterrows():
        print(f"  {row.Sym1}{row.Atom1} - {row.Sym2}{row.Atom2} : order {row.BondOrder}")
    print()

    return bond_table
