import os
import warnings

import numpy as np

from ._common import Struct
from .charges import _charge_color


def g16_charges_fchk(mol, ch, mode="atom", plot=True, atom_scale=0.35, bond_tol=1.30,
                     font_size=8, color_scale="RdBu", threshold=0.0):
    """Visualises atomic charges directly from a g16_fchk_read struct,
    without requiring the corresponding .log/.out file. Mirrors the
    interface of g16_charges exactly, so the two can be used
    interchangeably once mol/ch (i.e. data.mol/data.ch from
    g16_fchk_read) are available.

    Parameters
    ----------
    mol : Struct -- geometry struct from g16_fchk_read (data.mol);
        required fields: symbols, xyz, Natoms
    ch : Struct -- charge struct from g16_fchk_read (data.ch); required
        fields: charges, symbols, type, Natoms; optional: charges_H
        (.fchk files have no native H-summed data -- use g16_charges on
        the .log file for that)
    mode : 'atom' (default) | 'heavy'
    plot, atom_scale, bond_tol, font_size, color_scale, threshold : see
        g16_charges

    Returns
    -------
    ch_out : Struct -- symbols, charges, charges_H, sum_q, type, label,
        Natoms, filename (same layout as g16_charges).
    """
    mode = mode.lower()

    for field in ("symbols", "xyz", "Natoms"):
        if not hasattr(mol, field):
            raise ValueError(f'g16_charges_fchk: mol is missing field "{field}". '
                              "Use data.mol from g16_fchk_read.")
    for field in ("charges", "symbols", "type", "Natoms"):
        if not hasattr(ch, field):
            raise ValueError(f'g16_charges_fchk: ch is missing field "{field}". '
                              "Use data.ch from g16_fchk_read.")
    if mol.Natoms != ch.Natoms:
        raise ValueError(f"g16_charges_fchk: mol.Natoms ({mol.Natoms}) != ch.Natoms ({ch.Natoms}).")

    q_atom = np.asarray(ch.charges, dtype=float).ravel()
    q_heavy = None
    if getattr(ch, "charges_H", None) is not None and len(ch.charges_H) > 0:
        q_heavy = np.asarray(ch.charges_H, dtype=float).ravel()

    if mode == "atom":
        syms_use, xyz_use, q_use = mol.symbols, mol.xyz, q_atom
    elif mode == "heavy":
        if q_heavy is None:
            warnings.warn(
                "g16_charges_fchk: H-summed charges (charges_H) are not available in "
                ".fchk files. Falling back to per-atom charges. For H-summed charges, "
                "use g16_charges on the .log file."
            )
            syms_use, xyz_use, q_use = mol.symbols, mol.xyz, q_atom
        else:
            is_heavy = np.array([s != "H" for s in mol.symbols])
            syms_use = [s for s, h in zip(mol.symbols, is_heavy) if h]
            xyz_use = mol.xyz[is_heavy]
            q_use = q_heavy
            if len(q_use) != xyz_use.shape[0]:
                xyz_use = xyz_use[: len(q_use)]
    else:
        raise ValueError("g16_charges_fchk: mode must be 'atom' or 'heavy'.")

    src = getattr(mol, "filename", "") or getattr(ch, "filename", "")
    charge_type = ch.type
    found_label = getattr(ch, "label", "") or f"{charge_type} Charges (from .fchk)"

    ch_out = Struct(
        symbols=list(mol.symbols), charges=q_atom, charges_H=q_heavy,
        sum_q=float(q_atom.sum()), type=charge_type, label=found_label,
        Natoms=mol.Natoms, filename=src,
    )

    fname = os.path.splitext(os.path.basename(src))[0] if src else ""
    print(f"\n-- g16_charges_fchk ({charge_type}, {mode}): {fname} --")
    print(f"  Source : {found_label}")
    print(f"  {'Idx':>4}  {'Sym':<4}  {'q (e)':>8}")
    print(f"  {'-'*22}")
    for i, (s, q) in enumerate(zip(syms_use, q_use), 1):
        print(f"  {i:4d}  {s:<4}  {q:+8.4f}")
    print(f"  {'-'*22}")
    print(f"  Sum = {q_use.sum():+.5f} e\n")

    if not plot:
        return ch_out

    import matplotlib.pyplot as plt
    from .draw_molecule import g16_draw_molecule

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    g16_draw_molecule(mol, ax=ax, atom_scale=atom_scale, bond_tol=bond_tol,
                       show_labels=False, show_legend=True,
                       title=f"{fname} — {charge_type} charges ({mode})")

    q_max = np.max(np.abs(q_use)) if q_use.size else 1.0
    if q_max == 0:
        q_max = 1.0

    for i in range(len(q_use)):
        if abs(q_use[i]) < threshold:
            continue
        clr = _charge_color(q_use[i] / q_max) if color_scale.lower() == "rdbu" else (0.05, 0.05, 0.05)
        r_off = atom_scale * 0.8 + 0.3
        ax.text(xyz_use[i, 0], xyz_use[i, 1], xyz_use[i, 2] + r_off,
                f"{q_use[i]:+.3f}", fontsize=font_size, color=clr, fontweight="bold",
                ha="center", va="bottom")

    return ch_out
