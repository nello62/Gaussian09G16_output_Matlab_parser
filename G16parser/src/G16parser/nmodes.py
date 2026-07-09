import re

import numpy as np

from ._common import Struct, read_lines

_END_RE = re.compile(r"^\s*(-{20,}|Thermochemistry|Zero-point|Normal termination|Leave Link)")
_FREQ_RE = re.compile(r"^\s*Frequencies\s*--")
_REDMASS_RE = re.compile(r"^\s*Red\. masses\s*--")
_FRC_RE = re.compile(r"^\s*Frc consts\s*--")
_IR_RE = re.compile(r"^\s*IR Inten\s*--")
_RAMAN_RE = re.compile(r"^\s*Raman Activ\s*--")
_ATOM_HDR_RE = re.compile(r"^\s*Atom\s+AN\s+X")
_SYM_TOKEN_RE = re.compile(r"^[A-Za-z]")


def _parse_rhs(ln):
    idx = ln.find("--")
    return ln[idx + 2:] if idx != -1 else ln


def _floats(s):
    return [float(x) for x in s.split()]


def g16_nmodes(filename, section="last", modes=None, lines=None):
    """Extracts vibrational normal-mode displacement vectors from a
    Gaussian 16 .out/.log file.

    Parameters
    ----------
    filename : str
    section : 'last' (default) | 'first' — which "and normal coordinates:"
        section to read (opt+freq jobs can print more than one).
    modes : list[int], optional — 1-based mode indices to keep (default:
        all modes).
    lines : list[str], optional — pre-read file lines (see g16_read_all)

    Returns
    -------
    nm : Struct — freq, IR, Raman (np.ndarray, empty if absent), redmass,
        frcconst, symmetry (list[str]), disp (np.ndarray, shape
        (Natoms, 3, Nmodes)), Nmodes, Natoms, has_Raman, filename.
    """
    if lines is None:
        lines = read_lines(filename)
    n = len(lines)

    idx_sec = [i for i, ln in enumerate(lines) if "and normal coordinates:" in ln]
    if not idx_sec:
        raise ValueError(f'g16_nmodes: "normal coordinates" section not found in {filename}')

    sec_req = section.lower()
    if sec_req == "last":
        sec_start = idx_sec[-1]
    elif sec_req == "first":
        sec_start = idx_sec[0]
    else:
        raise ValueError("g16_nmodes: section must be 'first' or 'last'.")

    freqs, IRs, Ramans, redmass, frcconst = [], [], [], [], []
    symmetry = []
    disp_cols = []  # list of (Natoms, 3) arrays, one per mode column found
    natoms_det = 0

    k = sec_start + 1
    while k < n:
        ln = lines[k]

        if _END_RE.match(ln):
            break

        if _FREQ_RE.match(ln):
            if k > 0:
                ln_prev = lines[k - 1].strip()
                if "--" not in ln_prev and ln_prev:
                    for tok in ln_prev.split():
                        if _SYM_TOKEN_RE.match(tok):
                            symmetry.append(tok)
            freqs.extend(_floats(_parse_rhs(ln)))
            k += 1
            continue

        if _REDMASS_RE.match(ln):
            redmass.extend(_floats(_parse_rhs(ln)))
            k += 1
            continue

        if _FRC_RE.match(ln):
            frcconst.extend(_floats(_parse_rhs(ln)))
            k += 1
            continue

        if _IR_RE.match(ln):
            IRs.extend(_floats(_parse_rhs(ln)))
            k += 1
            continue

        if _RAMAN_RE.match(ln):
            Ramans.extend(_floats(_parse_rhs(ln)))
            k += 1
            continue

        if _ATOM_HDR_RE.match(ln):
            n_so_far = len(freqs)
            ncols = n_so_far if not disp_cols else n_so_far - len(disp_cols)
            if ncols < 1 or ncols > 3:
                ncols = 3

            k += 1
            atom_disp = []
            while k < n:
                tok2 = _floats(lines[k])
                if len(tok2) < 2 + 3 * ncols:
                    break
                atom_disp.append(tok2[2:2 + 3 * ncols])
                k += 1

            if not atom_disp:
                continue

            atom_disp = np.array(atom_disp)
            nat = atom_disp.shape[0]
            if natoms_det == 0:
                natoms_det = nat

            for ci in range(ncols):
                disp_cols.append(atom_disp[:, ci * 3:(ci + 1) * 3])
            continue

        k += 1

    nmodes = len(freqs)
    if nmodes == 0:
        raise ValueError(f"g16_nmodes: no modes read from {filename}")

    def fix_vec(v):
        v = np.array(v[:nmodes], dtype=float)
        if v.size < nmodes:
            v = np.concatenate([v, np.zeros(nmodes - v.size)])
        return v

    freqs = np.array(freqs, dtype=float)
    IRs = fix_vec(IRs)
    redmass = fix_vec(redmass)
    frcconst = fix_vec(frcconst)

    has_raman = len(Ramans) == nmodes
    Ramans = np.array(Ramans, dtype=float) if has_raman else np.array([])

    while len(symmetry) < nmodes:
        symmetry.append("?")
    symmetry = symmetry[:nmodes]

    disp_all = (np.stack(disp_cols, axis=2) if disp_cols
                else np.zeros((natoms_det, 3, 0)))

    if modes is not None:
        sel = [m for m in modes if 1 <= m <= nmodes]
        idx0 = [m - 1 for m in sel]
        freqs = freqs[idx0]
        IRs = IRs[idx0]
        if has_raman:
            Ramans = Ramans[idx0]
        redmass = redmass[idx0]
        frcconst = frcconst[idx0]
        symmetry = [symmetry[i] for i in idx0]
        if disp_all.shape[2] > 0:
            disp_all = disp_all[:, :, idx0]
        nmodes = len(freqs)

    nm = Struct(
        freq=freqs, IR=IRs, Raman=Ramans, redmass=redmass, frcconst=frcconst,
        symmetry=symmetry, disp=disp_all, Nmodes=nmodes, Natoms=natoms_det,
        has_Raman=has_raman, filename=filename,
    )
    print(f"g16_nmodes: {nmodes} modes, {natoms_det} atoms read from {filename}")
    return nm
