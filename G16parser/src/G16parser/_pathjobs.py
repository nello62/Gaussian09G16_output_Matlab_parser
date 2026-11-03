"""Shared helpers for multi-structure Gaussian 16 outputs (relaxed scans,
IRCs, --Link1-- chained jobs), used by g16_scan_read, g16_irc_read and
the multi-job check in g16_read_all.

Not part of the public API.
"""

import re

import numpy as np

from ._common import _SYMBOLS

HARTREE_TO_KCALMOL = 627.5095

SCF_RE = re.compile(r"SCF Done:\s+E\((\S+)\)\s*=\s*([-\d.]+)")
ORI_RE = re.compile(r"(Standard|Input) orientation\s*:", re.IGNORECASE)
_ROW6_RE = re.compile(r"^\s*(\d+)\s+(\d+)\s+(-?\d+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$")
_ROW5_RE = re.compile(r"^\s*(\d+)\s+(\d+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$")
_SYM_TO_Z = {s.upper(): i + 1 for i, s in enumerate(_SYMBOLS)}


def split_jobs(lines):
    """Returns a list of (start, end) 0-based line ranges, one per Gaussian
    job in the file. Jobs chained with --Link1-- each begin with an
    'Entering Link 1' line; the Freq step that Gaussian appends internally
    to an Opt+Freq route does NOT, so it stays in the same job as its Opt.
    """
    starts = [i for i, ln in enumerate(lines) if "Entering Link 1" in ln]
    if not starts:
        return [(0, len(lines))]
    starts[0] = 0
    ends = starts[1:] + [len(lines)]
    return list(zip(starts, ends))


def job_status(lines, a, b):
    """(completed, message) for the job in lines[a:b]."""
    for k in range(b - 1, a - 1, -1):
        ln = lines[k]
        if "Normal termination" in ln:
            return True, "Normal termination"
        if "Error termination" in ln:
            # the line(s) just before usually carry the actual reason
            msg = ""
            for j in range(k - 1, max(a, k - 6) - 1, -1):
                s = lines[j].strip()
                if s and not s.startswith(("Error", "Job cpu", "Elapsed", "File lengths")):
                    msg = s
                    break
            return False, "Error termination" + (f": {msg}" if msg else "")
    return False, "No termination line (job still running or file truncated)"


def orientation_label(lines, a, b):
    """'Standard orientation' if any such block exists in lines[a:b],
    otherwise 'Input orientation' (NoSymm jobs)."""
    for k in range(a, b):
        if "Standard orientation" in lines[k]:
            return "Standard orientation"
    return "Input orientation"


def read_orientation_block(lines, k):
    """Parses the 'Input/Standard orientation' block whose header is at
    line k. Returns (Z ndarray[int], xyz ndarray (N,3))."""
    Z, xyz = [], []
    for ln in lines[k + 5:]:
        s = ln.strip()
        if s and set(s) == {"-"}:
            break
        m = _ROW6_RE.match(ln)
        if m:
            Z.append(int(m.group(2)))
            xyz.append([float(m.group(4)), float(m.group(5)), float(m.group(6))])
    return np.array(Z, dtype=int), np.array(xyz, dtype=float).reshape(-1, 3)


def read_cartesian_block(lines, k):
    """Parses a 5-column 'Center / Atomic Number / X Y Z' block (as in the
    IRC 'CURRENT STRUCTURE' printout) starting the search at line k.
    Returns (Z, xyz) or (None, None) if no rows are found within 8 lines."""
    Z, xyz, started = [], [], False
    for j in range(k, len(lines)):
        ln = lines[j]
        m = _ROW5_RE.match(ln)
        if m:
            started = True
            Z.append(int(m.group(2)))
            xyz.append([float(m.group(3)), float(m.group(4)), float(m.group(5))])
            continue
        if started:
            break
        if j - k > 8:
            break
    if not Z:
        return None, None
    return np.array(Z, dtype=int), np.array(xyz, dtype=float).reshape(-1, 3)


def read_checkpoint_structure(lines, a, b):
    """Parses the 'Structure from the checkpoint file' block (lines like
    'C,0,x,y,z' or 'C,x,y,z') printed by Geom=Check jobs. Returns (Z, xyz)
    or (None, None)."""
    for k in range(a, b):
        if "Structure from the checkpoint file" in lines[k]:
            Z, xyz = [], []
            for ln in lines[k + 1:b]:
                parts = [p.strip() for p in ln.strip().split(",")]
                if len(parts) in (4, 5) and parts[0][:1].isalpha():
                    try:
                        x, y, z = (float(p) for p in parts[-3:])
                    except ValueError:
                        if Z:
                            break
                        continue
                    sym = re.sub(r"[^A-Za-z]", "", parts[0]).upper()
                    if sym not in _SYM_TO_Z:
                        break
                    Z.append(_SYM_TO_Z[sym])
                    xyz.append([x, y, z])
                elif Z:
                    break
            if Z:
                return np.array(Z, dtype=int), np.array(xyz, dtype=float)
    return None, None


def last_before(indices, k):
    """Largest element of the sorted list `indices` that is < k, or None."""
    best = None
    for i in indices:
        if i < k:
            best = i
        else:
            break
    return best


# ---------------------------------------------------------------------------
# internal coordinates (Angstrom, degrees) -- 0-based atom indices
# ---------------------------------------------------------------------------
def bond(xyz, i, j):
    return float(np.linalg.norm(xyz[i] - xyz[j]))


def angle(xyz, i, j, k):
    u, v = xyz[i] - xyz[j], xyz[k] - xyz[j]
    c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def dihedral(xyz, i, j, k, l):
    b0, b1, b2 = xyz[i] - xyz[j], xyz[k] - xyz[j], xyz[l] - xyz[k]
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def measure(xyz, ctype, atoms):
    """Value of a ModRedundant coordinate ('B','A','D') for 1-based atoms."""
    a = [n - 1 for n in atoms]
    if ctype == "B":
        return bond(xyz, *a)
    if ctype == "A":
        return angle(xyz, *a)
    if ctype == "D":
        return dihedral(xyz, *a)
    return float("nan")
