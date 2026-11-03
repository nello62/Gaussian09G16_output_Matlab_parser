import re
import warnings

import numpy as np

from ._common import Struct, read_lines, z_to_symbol, HARTREE_TO_KJMOL
from ._pathjobs import (HARTREE_TO_KCALMOL, SCF_RE, ORI_RE, split_jobs, job_status,
                        orientation_label, read_orientation_block, last_before, measure)

_SCANPT_RE = re.compile(r"on scan point\s+(\d+)\s+out of\s+(\d+)")
_STEP_RE = re.compile(r"Step number\s+(\d+)\s+out of a maximum of\s+\d+\s+on scan point\s+(\d+)")
_MODRED_RE = re.compile(r"^\s*([BADL])\s+((?:\d+\s+){2,4})([SF])\b(?:\s+(\d+)\s+([-\d.]+))?",
                        re.IGNORECASE)
_UNITS = {"B": "Angstrom", "A": "deg", "D": "deg", "L": "deg"}


def g16_scan_read(filename, job="auto", plot=False, ax=None, lines=None):
    """Reads a relaxed potential-energy-surface scan (Opt=ModRedundant with an 'S' coordinate) point by point.

    One entry is returned per scan point: the converged constrained
    geometry, its SCF energy, the scanned coordinate(s) measured on that
    geometry, and the relative energy. Works on files where the scan is
    one job of a --Link1-- chain (the scan job is located automatically).

    Parameters
    ----------
    filename : str
    job : 'auto' (default) | int
        Which job of a --Link1-- chain to read (1-based). 'auto' takes the
        first job containing a relaxed scan.
    plot : bool
        If True, plots dE (kcal/mol) versus the scanned coordinate (1-D
        scan) or versus the point index (multi-coordinate scan).
    ax : matplotlib Axes, optional
    lines : list[str], optional
        Pre-read file lines (see g16_read_all).

    Returns
    -------
    S : Struct with fields
        coords        list[Struct]  scanned coordinates: type ('B','A','D'),
                      atoms (1-based tuple), label (e.g. 'D(4,1,3,9)'),
                      nsteps, stepsize, units
        frozen        list[str]     frozen ('F') coordinates, as labels
        value         ndarray (Npts, Ncoord)  scanned coordinate(s) measured
                      on each geometry (dihedrals unwrapped along the scan)
        E             ndarray (Npts,)  SCF energy (Hartree)
        dE_kcal, dE_kJ  ndarray (Npts,)  energy relative to the lowest point
        i_min         int           0-based index of the lowest point
        point         ndarray (Npts,)  scan-point number (1-based)
        converged     ndarray (Npts,)  bool, 'Stationary point found'
        opt_steps     ndarray (Npts,)  optimisation steps used per point
        n_points_expected int       total points announced by Gaussian
        mol           list[Struct]  per-point structures, compatible with
                      g16_get_bond_length / g16_get_angle / g16_get_dihedral /
                      g16_draw_molecule
        xyz           ndarray (Npts, Natoms, 3)
        symbols, Z, Natoms, orientation, method
        completed     bool          job ended with Normal termination
        status        str
        job           int           1-based job index used
        filename      str

    Notes
    -----
    Gaussian prints the dipole moment and atomic charges of a relaxed
    scan only for the first and last geometry, so they are not returned
    per point; run single points on S.mol[k] if they are needed.
    """
    if lines is None:
        lines = read_lines(filename)

    jobs = split_jobs(lines)
    if job == "auto":
        cand = [n for n, (a, b) in enumerate(jobs)
                if any(_SCANPT_RE.search(ln) for ln in lines[a:b])]
        if not cand:
            raise ValueError(f"g16_scan_read: no relaxed scan ('on scan point') found in {filename}")
        jidx = cand[0]
    else:
        jidx = int(job) - 1
        if jidx < 0 or jidx >= len(jobs):
            raise ValueError(f"g16_scan_read: job {job} out of range [1, {len(jobs)}].")
    a, b = jobs[jidx]

    # ---- scanned / frozen coordinates from the ModRedundant echo
    coords, frozen = [], []
    for k in range(a, b):
        if "ModRedundant input section has been read" in lines[k]:
            for ln in lines[k + 1:b]:
                m = _MODRED_RE.match(ln)
                if not m:
                    break
                ctype = m.group(1).upper()
                atoms = tuple(int(t) for t in m.group(2).split())
                label = f"{ctype}({','.join(str(t) for t in atoms)})"
                if m.group(3).upper() == "S":
                    coords.append(Struct(type=ctype, atoms=atoms, label=label,
                                         nsteps=int(m.group(4)) if m.group(4) else None,
                                         stepsize=float(m.group(5)) if m.group(5) else None,
                                         units=_UNITS.get(ctype, "")))
                else:
                    frozen.append(label)
            break

    # ---- markers
    scf_idx, scf_val, method = [], [], ""
    ori_idx, scan_marks, steps = [], [], {}
    stat_idx = []
    ori_label = orientation_label(lines, a, b)
    n_expected = 0
    for k in range(a, b):
        ln = lines[k]
        m = SCF_RE.search(ln)
        if m:
            scf_idx.append(k)
            scf_val.append(float(m.group(2)))
            method = m.group(1)
            continue
        if ori_label in ln and ORI_RE.search(ln):
            ori_idx.append(k)
            continue
        m = _SCANPT_RE.search(ln)
        if m:
            scan_marks.append((k, int(m.group(1))))
            n_expected = int(m.group(2))
            ms = _STEP_RE.search(ln)
            if ms:
                p = int(ms.group(2))
                steps[p] = max(steps.get(p, 0), int(ms.group(1)))
            continue
        if "Stationary point found" in ln:
            stat_idx.append(k)

    if not scan_marks:
        raise ValueError(f"g16_scan_read: job {jidx + 1} of {filename} is not a relaxed scan.")

    mark_lines = [k for k, _ in scan_marks]

    def point_at(k):
        j = last_before(mark_lines, k)
        return None if j is None else dict(scan_marks)[j]

    ends = []                                   # (line index, converged)
    seen = set()
    for k in stat_idx:
        p = point_at(k)
        if p is not None and p not in seen:
            ends.append((k, p, True))
            seen.add(p)
    last_p = scan_marks[-1][1]
    if last_p not in seen:                      # job stopped inside this point
        ends.append((b, last_p, False))

    E, pts, conv, xyz_list, Zref = [], [], [], [], None
    for k, p, ok in ends:
        i_scf = last_before(scf_idx, k)
        i_ori = last_before(ori_idx, k)
        if i_scf is None or i_ori is None:
            continue
        Z, xyz = read_orientation_block(lines, i_ori)
        Zref = Z if Zref is None else Zref
        E.append(scf_val[scf_idx.index(i_scf)])
        pts.append(p)
        conv.append(ok)
        xyz_list.append(xyz)

    if not E:
        raise ValueError(f"g16_scan_read: no scan point could be read from {filename}")

    E = np.array(E)
    xyz_arr = np.stack(xyz_list)
    symbols = [z_to_symbol(z) for z in Zref]
    value = np.full((len(E), max(len(coords), 1)), np.nan)
    for c, co in enumerate(coords):
        v = np.array([measure(x, co.type, co.atoms) for x in xyz_arr])
        if co.type in ("D", "L"):
            v = np.degrees(np.unwrap(np.radians(v)))
        value[:, c] = v
    if not coords:
        value = value[:, :0]

    i_min = int(np.argmin(E))
    dE_h = E - E[i_min]
    completed, status = job_status(lines, a, b)
    if not all(conv):
        warnings.warn(f"g16_scan_read: {sum(not c for c in conv)} scan point(s) did not converge "
                      f"in {filename} ({status}).")

    mol = [Struct(symbols=symbols, xyz=x, Z=Zref.copy(), Natoms=len(Zref), step=n + 1,
                  n_steps=len(E), orientation=ori_label, filename=filename)
           for n, x in enumerate(xyz_arr)]

    S = Struct(
        coords=coords, frozen=frozen, value=value, E=E,
        dE_kcal=dE_h * HARTREE_TO_KCALMOL, dE_kJ=dE_h * HARTREE_TO_KJMOL, i_min=i_min,
        point=np.array(pts), converged=np.array(conv),
        opt_steps=np.array([steps.get(p, 0) for p in pts]),
        n_points_expected=n_expected, mol=mol, xyz=xyz_arr,
        symbols=symbols, Z=Zref, Natoms=len(Zref), orientation=ori_label, method=method,
        completed=completed, status=status, job=jidx + 1, filename=filename,
    )

    lab = ", ".join(c.label for c in coords) if coords else "(no S coordinate echoed)"
    print(f"\n-- g16_scan_read: {filename} (job {jidx + 1}) --")
    print(f"  Scanned: {lab}   points read: {len(E)}/{n_expected}   {status}")
    print(f"  E(min) = {E[i_min]:.8f} Ha at point {pts[i_min]};  dE(max) = {S.dE_kcal.max():.3f} kcal/mol\n")

    if plot:
        _plot_scan(S, ax)
    return S


def _plot_scan(S, ax=None):
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(6.5, 4.0))
    if len(S.coords) == 1:
        x = S.value[:, 0]
        ax.set_xlabel(f"{S.coords[0].label} ({S.coords[0].units})")
    else:
        x = S.point
        ax.set_xlabel("Scan point")
    ax.plot(x, S.dE_kcal, "o-", ms=4)
    bad = ~S.converged
    if bad.any():
        ax.plot(x[bad], S.dE_kcal[bad], "x", color="red", ms=8, label="not converged")
        ax.legend()
    ax.set_ylabel("ΔE (kcal/mol)")
    ax.set_title(f"Relaxed scan — {S.method}")
    ax.grid(alpha=0.3)
    return ax
