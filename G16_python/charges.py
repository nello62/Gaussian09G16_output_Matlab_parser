import math
import os
import re
import warnings

import numpy as np

from ._common import Struct, read_lines, DEBYE_TO_AU
from .dipole_polar import g16_dipole_polar

_ATOM_ROW_RE = re.compile(r"^\s*\d+\s+([A-Za-z]+)\s+(-?[\d.]+)")


def _parse_block(lines, k_start, n):
    syms, qs = [], []
    k2 = k_start + 2  # skip the "    1" index row
    while k2 < n:
        ln2 = lines[k2].strip()
        if not ln2:
            break
        low = ln2.lower()
        if "sum of" in low or "charges" in low:
            break
        m = _ATOM_ROW_RE.match(lines[k2])
        if m:
            syms.append(m.group(1))
            qs.append(float(m.group(2)))
        k2 += 1
    return syms, np.array(qs)


def _charge_color(t):
    t = max(-1.0, min(1.0, t))
    if t >= 0:
        c = (1, 1 - t, 1 - t)
    else:
        c = (1 + t, 1 + t, 1)
    return tuple(x * 0.82 for x in c)


def _charge_centroid(xyz, charges, mask, label):
    if not np.any(mask):
        warnings.warn(f"No atoms with a {label} partial charge; using the unweighted atom centroid.")
        return xyz.mean(axis=0)
    w = np.abs(charges[mask])
    return (xyz[mask] * w[:, None]).sum(axis=0) / w.sum()


def _dipole_origin(spec, xyz, charges):
    if isinstance(spec, (tuple, list, np.ndarray)) and len(spec) == 3:
        return np.array(spec, dtype=float)
    spec_low = spec.lower()
    if spec_low == "negcharge":
        return _charge_centroid(xyz, charges, charges < 0, "negative")
    if spec_low == "poscharge":
        return _charge_centroid(xyz, charges, charges > 0, "positive")
    if spec_low == "centroid":
        return xyz.mean(axis=0)
    warnings.warn(f"Unknown DipoleOrigin = \"{spec}\"; using the unweighted atom centroid.")
    return xyz.mean(axis=0)


def _origin_label(spec):
    if isinstance(spec, (tuple, list, np.ndarray)):
        return str(np.round(np.asarray(spec), 3).tolist())
    return spec


def g16_charges(filename, type="Mulliken", mode="atom", plot=True,
                atom_scale=0.35, bond_tol=1.30, font_size=8, color_scale="RdBu",
                threshold=0.0, show_dipole=False, dipole_origin="negcharge",
                dipole_scale=1.0, dipole_color=(0, 0.6, 0), dipole_line_width=2.5,
                show_dipole_label=True, dipole_font_size=11, dipole_units="Debye",
                lines=None):
    """Extracts Mulliken or APT atomic charges from a Gaussian 16 output
    file, and optionally renders them on the 3D structure with a dipole
    moment arrow overlay.

    See G16_charges.m for the full parameter/field documentation — this
    is a line-by-line Python port with the same defaults and behaviour.

    Returns
    -------
    ch : Struct — symbols, charges, charges_H, sum_q, type, label, Natoms,
        filename, dipole, dipole_origin, dipole_Debye, dipole_au.
    """
    mode = mode.lower()
    dipole_units = dipole_units.lower()
    if dipole_units not in ("debye", "au"):
        warnings.warn(f"Unknown DipoleUnits = \"{dipole_units}\"; using 'Debye'.")
        dipole_units = "debye"

    if lines is None:
        lines = read_lines(filename)
    n = len(lines)

    atom_starts, heavy_starts, found_labels = [], [], []
    for i, ln in enumerate(lines):
        low = ln.lower()
        if type.lower() not in low:
            continue
        if "charges" not in low:
            continue
        trimmed = ln.strip()
        if not trimmed or not trimmed.endswith(":"):
            continue
        if "sum of" in low:
            continue
        if "hydrogen" in low or "summed" in low:
            heavy_starts.append(i)
        else:
            atom_starts.append(i)
            found_labels.append(trimmed)

    if not atom_starts:
        candidates = [f"  line {i+1}: {lines[i].strip()}"
                      for i in range(min(n, 300)) if type.lower() in lines[i].lower()]
        if not candidates:
            raise ValueError(f'g16_charges: no "{type}" charges found in {filename}')
        raise ValueError(
            f'g16_charges: charge header not found in {filename}\n"{type}" appears in:\n'
            + "\n".join(candidates)
        )

    syms_atom, q_atom = _parse_block(lines, atom_starts[-1], n)
    if q_atom.size == 0:
        raise ValueError(f"g16_charges: charge header found but no atom data read from {filename}")
    found_label = found_labels[-1]

    syms_heavy, q_heavy = [], np.array([])
    if heavy_starts:
        syms_heavy, q_heavy = _parse_block(lines, heavy_starts[-1], n)

    if mode == "atom":
        syms_use, q_use = syms_atom, q_atom
    elif mode == "heavy":
        if q_heavy.size == 0:
            warnings.warn("g16_charges: H-summed charges not found; falling back to per-atom.")
            syms_use, q_use = syms_atom, q_atom
        else:
            syms_use, q_use = syms_heavy, q_heavy
    else:
        raise ValueError("g16_charges: mode must be 'atom' or 'heavy'.")

    ch = Struct(
        symbols=syms_atom, charges=q_atom, charges_H=q_heavy, sum_q=float(q_atom.sum()),
        type=type, label=found_label, Natoms=len(q_atom), filename=filename,
        dipole=None, dipole_origin=None, dipole_Debye=None, dipole_au=None,
    )

    mol = None
    if show_dipole:
        from .structure import g16_structure
        mu = None
        try:
            dp = g16_dipole_polar(filename, units="debye", lines=lines)
            mu = np.array([dp.mu_x, dp.mu_y, dp.mu_z])
        except Exception as exc:
            warnings.warn(f"Could not read the dipole moment ({exc}); ShowDipole will be ignored.")

        if mu is None or np.linalg.norm(mu) < np.finfo(float).eps:
            if mu is None:
                warnings.warn("Dipole moment field not recognised in g16_dipole_polar output; ShowDipole will be ignored.")
        else:
            mol = g16_structure(filename, lines=lines)
            ch.dipole = mu
            ch.dipole_origin = _dipole_origin(dipole_origin, mol.xyz, q_atom)
            ch.dipole_Debye = float(np.linalg.norm(mu))
            ch.dipole_au = ch.dipole_Debye * DEBYE_TO_AU

    dip_value_disp = dip_unit_symbol = None
    if ch.dipole is not None:
        if dipole_units == "au":
            dip_value_disp, dip_unit_symbol = ch.dipole_au, "a.u."
        else:
            dip_value_disp, dip_unit_symbol = ch.dipole_Debye, "D"

    print(f"\n── g16_charges ({type}, {mode}): {filename} ──")
    print(f'  Header found : "{found_label}"')
    print(f"  {'Idx':>4}  {'Sym':<4}  {'q (e)':>8}")
    print(f"  {'-'*22}")
    for i, (s, q) in enumerate(zip(syms_use, q_use), 1):
        print(f"  {i:4d}  {s:<4}  {q:+8.4f}")
    print(f"  {'-'*22}")
    print(f"  Sum = {q_use.sum():+.5f} e")
    if ch.dipole is not None:
        print(f"  |mu| = {dip_value_disp:.3f} {dip_unit_symbol}  (anchor: {_origin_label(dipole_origin)})")
    print()

    if plot:
        import matplotlib.pyplot as plt
        from .draw_molecule import g16_draw_molecule

        if mol is None:
            from .structure import g16_structure
            mol = g16_structure(filename, lines=lines)
        fname = os.path.splitext(os.path.basename(filename))[0]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        g16_draw_molecule(mol, ax=ax, atom_scale=atom_scale, bond_tol=bond_tol,
                           show_labels=False, show_legend=True,
                           title=f"{fname} — {type} charges ({mode})")

        if mode == "atom":
            xyz_use = mol.xyz
        else:
            is_heavy = np.array([s != "H" for s in mol.symbols])
            xyz_use = mol.xyz[is_heavy]
            if xyz_use.shape[0] != len(q_use):
                xyz_use = mol.xyz[: len(q_use)]

        q_max = np.max(np.abs(q_use)) if q_use.size else 1.0
        if q_max == 0:
            q_max = 1.0

        for i in range(len(q_use)):
            if abs(q_use[i]) < threshold:
                continue
            clr = _charge_color(q_use[i] / q_max) if color_scale.lower() == "rdbu" else (0, 0, 0)
            r_off = atom_scale * 0.8 + 0.3
            ax.text(xyz_use[i, 0], xyz_use[i, 1], xyz_use[i, 2] + r_off,
                    f"{q_use[i]:+.3f}", fontsize=font_size, color=clr,
                    fontweight="bold", ha="center", va="bottom")

        if show_dipole and ch.dipole is not None:
            mu = ch.dipole
            origin = ch.dipole_origin
            vec = mu * dipole_scale  # cosmetic length only (already Debye-scaled)
            ax.quiver(origin[0], origin[1], origin[2], vec[0], vec[1], vec[2],
                      color=dipole_color, linewidth=dipole_line_width)
            if show_dipole_label:
                tip = origin + vec
                ax.text(tip[0], tip[1], tip[2], f"  mu = {dip_value_disp:.2f} {dip_unit_symbol}",
                        color=dipole_color, fontsize=dipole_font_size, fontweight="bold")

        plt.show()

    return ch
