import math
import re

from ._common import Struct, read_lines, HARTREE_TO_KJMOL

_SCF_RE = re.compile(r"SCF Done:\s+E\((\S+)\)\s*=\s*(-?[\d.]+)")
_TEMP_P_RE = re.compile(r"Temperature\s+([\d.]+)\s+Kelvin\.\s+Pressure\s+([\d.]+)")
_ZPE_KJ_RE = re.compile(r"Zero-point vibrational energy\s+([\d.]+)\s+\(Joules")
_ZPE_CORR_RE = re.compile(r"Zero-point correction=\s+(-?[\d.]+)")
_U_CORR_RE = re.compile(r"Thermal correction to Energy=\s+(-?[\d.]+)")
_H_CORR_RE = re.compile(r"Thermal correction to Enthalpy=\s+(-?[\d.]+)")
_G_CORR_RE = re.compile(r"Thermal correction to Gibbs Free Energy=\s+(-?[\d.]+)")
_E0_RE = re.compile(r"Sum of electronic and zero-point Energies=\s+(-?[\d.]+)")
_U_RE = re.compile(r"Sum of electronic and thermal Energies=\s+(-?[\d.]+)")
_H_RE = re.compile(r"Sum of electronic and thermal Enthalpies=\s+(-?[\d.]+)")
_G_RE = re.compile(r"Sum of electronic and thermal Free Energies=\s+(-?[\d.]+)")


def g16_energy(filename, step="last", lines=None):
    """Extracts SCF energy and thermochemistry from a Gaussian 16 .out/.log file.

    Parameters
    ----------
    filename : str
    step : 'last' (default) | 'first' | int (1-based) — which SCF Done
        occurrence to report.
    lines : list[str], optional
        Pre-read file lines, to skip re-reading the file when it has
        already been read elsewhere (see g16_read_all).

    Returns
    -------
    en : Struct — see G16_energy.m for the full field list (SCF, method,
        ZPE_corr, U_corr, H_corr, G_corr, E0, U, H, G, ZPE_kJ, T, P,
        has_thermo, S_JmolK, SCF_kJ, G_kJ, H_kJ, filename). For opt-only
        jobs (no freq) has_thermo is False and the thermo fields are NaN.
    """
    if lines is None:
        lines = read_lines(filename)

    scf_vals, scf_methods = [], []
    for ln in lines:
        m = _SCF_RE.search(ln)
        if m:
            scf_methods.append(m.group(1))
            scf_vals.append(float(m.group(2)))

    if not scf_vals:
        raise ValueError(f'g16_energy: no "SCF Done" line found in {filename}')

    n_scf = len(scf_vals)
    if isinstance(step, str):
        if step.lower() == "last":
            si = n_scf
        elif step.lower() == "first":
            si = 1
        else:
            raise ValueError("g16_energy: step must be 'first', 'last' or an integer.")
    else:
        si = round(step)
        if si < 1 or si > n_scf:
            raise ValueError(f"g16_energy: step {si} out of range [1,{n_scf}].")

    SCF = scf_vals[si - 1]
    method = scf_methods[si - 1]

    ZPE_corr = U_corr = H_corr = G_corr = math.nan
    E0 = U_tot = H_tot = G_tot = math.nan
    ZPE_kJ = T = P_atm = math.nan

    for ln in lines:
        m = _TEMP_P_RE.search(ln)
        if m:
            T, P_atm = float(m.group(1)), float(m.group(2))
            continue
        m = _ZPE_KJ_RE.search(ln)
        if m:
            ZPE_kJ = float(m.group(1)) / 1000  # J/mol -> kJ/mol
            continue
        m = _ZPE_CORR_RE.search(ln)
        if m:
            ZPE_corr = float(m.group(1)); continue
        m = _U_CORR_RE.search(ln)
        if m:
            U_corr = float(m.group(1)); continue
        m = _H_CORR_RE.search(ln)
        if m:
            H_corr = float(m.group(1)); continue
        m = _G_CORR_RE.search(ln)
        if m:
            G_corr = float(m.group(1)); continue
        m = _E0_RE.search(ln)
        if m:
            E0 = float(m.group(1)); continue
        m = _U_RE.search(ln)
        if m:
            U_tot = float(m.group(1)); continue
        m = _H_RE.search(ln)
        if m:
            H_tot = float(m.group(1)); continue
        m = _G_RE.search(ln)
        if m:
            G_tot = float(m.group(1)); continue

    has_thermo = not math.isnan(ZPE_corr)
    if has_thermo and not math.isnan(T) and T > 0:
        S_JmolK = (H_tot - G_tot) / T * HARTREE_TO_KJMOL * 1000
    else:
        S_JmolK = math.nan

    en = Struct(
        SCF=SCF, method=method,
        ZPE_corr=ZPE_corr, U_corr=U_corr, H_corr=H_corr, G_corr=G_corr,
        E0=E0, U=U_tot, H=H_tot, G=G_tot,
        ZPE_kJ=ZPE_kJ, T=T, P=P_atm, has_thermo=has_thermo,
        S_JmolK=S_JmolK,
        SCF_kJ=SCF * HARTREE_TO_KJMOL,
        G_kJ=(G_tot * HARTREE_TO_KJMOL if has_thermo else math.nan),
        H_kJ=(H_tot * HARTREE_TO_KJMOL if has_thermo else math.nan),
        filename=filename,
    )

    print(f"\n── g16_energy: {filename} ──")
    print(f"  Method : {en.method}")
    print(f"  SCF    : {en.SCF:+.8f}  Ha")
    if en.has_thermo:
        print(f"  ZPE    : {en.ZPE_corr:+.8f}  Ha   ({en.ZPE_kJ:.2f} kJ/mol)")
        print(f"  E0+ZPE : {en.E0:+.8f}  Ha")
        print(f"  H      : {en.H:+.8f}  Ha")
        print(f"  G      : {en.G:+.8f}  Ha")
        print(f"  S      :  {en.S_JmolK:.4f}  J/(mol*K)")
        print(f"  T = {en.T:.2f} K,  P = {en.P:.5f} atm")
    print()

    return en
