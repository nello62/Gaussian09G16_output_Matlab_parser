import re
import warnings

import numpy as np

from ._common import Struct, read_lines, z_to_symbol, HARTREE_TO_KJMOL
from ._pathjobs import (HARTREE_TO_KCALMOL, SCF_RE, ORI_RE, split_jobs, job_status,
                        read_orientation_block, read_cartesian_block,
                        read_checkpoint_structure, last_before)

_PT_RE = re.compile(r"Point Number:\s+(\d+)\s+Path Number:\s+(\d+)")
_DIR_RE = re.compile(r"Point Number\s+(\d+)\s+in\s+(FORWARD|REVERSE)\s+path direction", re.IGNORECASE)
_NET_RE = re.compile(r"NET REACTION COORDINATE UP TO THIS POINT\s*=\s*([-\d.]+)")
_CORR_RE = re.compile(r"Corrected End Point Energy\s*=\s*([-\d.]+)")
_CHK_E_RE = re.compile(r"Energy From Chk\s*=\s*([-\d.]+)")
_MAXPTS_RE = re.compile(r"Maximum points per path\s*=\s*(\d+)")
_STEP_RE = re.compile(r"Step size\s*=\s*([\d.]+)\s*(\S+)")
_SCHEME_RE = re.compile(r"Integration scheme\s*=\s*(\S+)")


def g16_irc_read(filename, job="auto", plot=False, ax=None, lines=None):
    """Reads an IRC (intrinsic reaction coordinate) calculation point by point, ordered reverse -> TS -> forward.

    Works on a stand-alone IRC job, on an IRC that is one job of a
    --Link1-- chain (e.g. TS optimisation + Freq + IRC), and on an IRC
    that stopped with an error: every point Gaussian completed is returned.

    Parameters
    ----------
    filename : str
    job : 'auto' (default) | int
        Which job of a --Link1-- chain to read (1-based). 'auto' takes the
        first job containing an IRC.
    plot : bool
        If True, plots dE (kcal/mol, relative to the TS) versus the
        reaction coordinate.
    ax : matplotlib Axes, optional
    lines : list[str], optional
        Pre-read file lines (see g16_read_all).

    Returns
    -------
    I : Struct with fields
        rxcoord     ndarray (Npts,)  reaction coordinate (sqrt(amu)*bohr),
                    negative on the REVERSE branch, 0 at the TS
        E           ndarray (Npts,)  energy (Hartree), NaN if not printed
        dE_kcal, dE_kJ  ndarray       energy relative to the TS
        E_source    list[str]        'chk' (TS energy read from checkpoint),
                    'scf' (SCF at the point) or 'hpc_corrected' (HPC
                    corrector end-point energy, DWI-interpolated)
        direction   list[str]        'REVERSE' | 'TS' | 'FORWARD'
        point       ndarray (Npts,)  Gaussian point number on its branch
        i_ts        int              0-based index of the TS
        n_forward, n_reverse  int    points completed on each branch
        mol         list[Struct]     per-point structures (g16_get_* / draw compatible)
        xyz         ndarray (Npts, Natoms, 3)
        symbols, Z, Natoms
        max_points  int | None       MaxPoints per path requested
        step_size   str              e.g. '0.100 bohr'
        scheme      str              integration scheme (HPC, LQA, ...)
        completed   bool             job ended with Normal termination
        status      str
        job         int              1-based job index used
        filename    str
    """
    if lines is None:
        lines = read_lines(filename)

    jobs = split_jobs(lines)
    if job == "auto":
        cand = [n for n, (a, b) in enumerate(jobs)
                if any("IRC-IRC-IRC" in ln for ln in lines[a:b])]
        if not cand:
            raise ValueError(f"g16_irc_read: no IRC calculation found in {filename}")
        jidx = cand[0]
    else:
        jidx = int(job) - 1
        if jidx < 0 or jidx >= len(jobs):
            raise ValueError(f"g16_irc_read: job {job} out of range [1, {len(jobs)}].")
    a, b = jobs[jidx]

    scf_idx, scf_val, ori_idx = [], [], []
    corr_idx, corr_val = [], []
    dir_marks = []                     # (line, point, direction)
    pt_marks = []                      # (line, point, path)
    max_points, step_size, scheme, E_chk = None, "", "", None
    for k in range(a, b):
        ln = lines[k]
        m = SCF_RE.search(ln)
        if m:
            scf_idx.append(k)
            scf_val.append(float(m.group(2)))
            continue
        if ORI_RE.search(ln):
            ori_idx.append(k)
            continue
        m = _CORR_RE.search(ln)
        if m:
            corr_idx.append(k)
            corr_val.append(float(m.group(1)))
            continue
        m = _PT_RE.search(ln)
        if m:
            pt_marks.append((k, int(m.group(1)), int(m.group(2))))
            continue
        m = _DIR_RE.search(ln)
        if m:
            dir_marks.append((k, int(m.group(1)), m.group(2).upper()))
            continue
        m = _CHK_E_RE.search(ln)
        if m and E_chk is None:
            E_chk = float(m.group(1))
            continue
        if max_points is None:
            m = _MAXPTS_RE.search(ln)
            if m:
                max_points = int(m.group(1))
                continue
        if not step_size:
            m = _STEP_RE.search(ln)
            if m and "Supplied" not in ln:
                step_size = f"{m.group(1)} {m.group(2)}"
                continue
        if not scheme:
            m = _SCHEME_RE.search(ln)
            if m:
                scheme = m.group(1)

    if not pt_marks:
        raise ValueError(f"g16_irc_read: no IRC point found in job {jidx + 1} of {filename}")

    dir_lines = [d[0] for d in dir_marks]
    rows = []                          # dict per point
    seen = set()
    prev_line = a
    for k, p, path in pt_marks:
        if p == 0:
            if any(r["direction"] == "TS" for r in rows):
                prev_line = k
                continue
            # TS geometry: last orientation block before point 0, else the
            # checkpoint structure echoed by Geom=Check
            i_ori = last_before(ori_idx, k)
            if i_ori is not None:
                Z, xyz = read_orientation_block(lines, i_ori)
            else:
                Z, xyz = read_checkpoint_structure(lines, a, k)
            if E_chk is not None:
                E, src = E_chk, "chk"
            else:
                i_scf = last_before(scf_idx, k)
                E, src = (scf_val[scf_idx.index(i_scf)], "scf") if i_scf is not None else (np.nan, "none")
            rows.append(dict(direction="TS", point=0, rx=0.0, E=E, src=src, Z=Z, xyz=xyz))
            prev_line = k
            continue

        j = last_before(dir_lines, k)
        direction = dir_marks[dir_lines.index(j)][2] if j is not None else ("FORWARD" if path == 1 else "REVERSE")
        if (direction, p) in seen:
            prev_line = k
            continue
        seen.add((direction, p))

        Z, xyz = read_cartesian_block(lines, k + 1)
        if Z is None:                  # fallback: last orientation block
            i_ori = last_before(ori_idx, k)
            Z, xyz = read_orientation_block(lines, i_ori) if i_ori is not None else (None, None)

        rx = np.nan
        for ln in lines[k + 1:min(b, k + 40)]:
            m = _NET_RE.search(ln)
            if m:
                rx = float(m.group(1))
                break

        i_corr = last_before(corr_idx, k)
        if i_corr is not None and i_corr > prev_line:
            E, src = corr_val[corr_idx.index(i_corr)], "hpc_corrected"
        else:
            i_scf = last_before(scf_idx, k)
            E, src = (scf_val[scf_idx.index(i_scf)], "scf") if i_scf is not None and i_scf > prev_line \
                else (np.nan, "none")

        sign = -1.0 if direction == "REVERSE" else 1.0
        rows.append(dict(direction=direction, point=p, rx=sign * abs(rx), E=E, src=src, Z=Z, xyz=xyz))
        prev_line = k

    rows = [r for r in rows if r["Z"] is not None]
    rev = sorted([r for r in rows if r["direction"] == "REVERSE"], key=lambda r: -r["point"])
    ts = [r for r in rows if r["direction"] == "TS"]
    fwd = sorted([r for r in rows if r["direction"] == "FORWARD"], key=lambda r: r["point"])
    ordered = rev + ts + fwd
    if not ts:
        warnings.warn(f"g16_irc_read: TS (point 0) geometry not found in {filename}; "
                      "energies are given relative to the highest point instead.")

    Zref = ordered[0]["Z"]
    symbols = [z_to_symbol(z) for z in Zref]
    E = np.array([r["E"] for r in ordered], dtype=float)
    i_ts = len(rev) if ts else int(np.nanargmax(E))
    dE_h = E - E[i_ts]
    xyz_arr = np.stack([r["xyz"] for r in ordered])
    completed, status = job_status(lines, a, b)

    mol = [Struct(symbols=symbols, xyz=r["xyz"], Z=Zref.copy(), Natoms=len(Zref), step=n + 1,
                  n_steps=len(ordered), orientation="IRC path", filename=filename)
           for n, r in enumerate(ordered)]

    I = Struct(
        rxcoord=np.array([r["rx"] for r in ordered]), E=E,
        dE_kcal=dE_h * HARTREE_TO_KCALMOL, dE_kJ=dE_h * HARTREE_TO_KJMOL,
        E_source=[r["src"] for r in ordered], direction=[r["direction"] for r in ordered],
        point=np.array([r["point"] for r in ordered]), i_ts=i_ts,
        n_forward=len(fwd), n_reverse=len(rev), mol=mol, xyz=xyz_arr,
        symbols=symbols, Z=Zref, Natoms=len(Zref),
        max_points=max_points, step_size=step_size, scheme=scheme,
        completed=completed, status=status, job=jidx + 1, filename=filename,
    )

    print(f"\n-- g16_irc_read: {filename} (job {jidx + 1}) --")
    print(f"  Points: {len(rev)} reverse + {'TS' if ts else 'no TS'} + {len(fwd)} forward"
          f"   (MaxPoints={max_points}, step {step_size}, {scheme})")
    print(f"  {status}")
    if not completed:
        print("  Note: the IRC stopped before the end; all completed points are returned.")
    print()

    if plot:
        _plot_irc(I, ax)
    return I


def _plot_irc(I, ax=None):
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.plot(I.rxcoord, I.dE_kcal, "o-", ms=3)
    ax.plot(I.rxcoord[I.i_ts], I.dE_kcal[I.i_ts], "*", ms=12, color="crimson", label="TS")
    ax.axvline(0, color="0.7", lw=0.8, ls=":")
    ax.set_xlabel("Reaction coordinate (amu$^{1/2}$ bohr)  ← reverse | forward →")
    ax.set_ylabel("ΔE relative to TS (kcal/mol)")
    ax.set_title("IRC" + ("" if I.completed else " (incomplete)"))
    ax.legend()
    ax.grid(alpha=0.3)
    return ax
