import math
import re

import numpy as np

from ._common import Struct, read_lines, fortran_to_float, AU_TO_DEBYE

_COMPACT_DIPOLE_RE = re.compile(
    r"X=\s*(-?[\d.]+)\s+Y=\s*(-?[\d.]+)\s+Z=\s*(-?[\d.]+)\s+Tot=\s*([\d.]+)"
)
_COMPACT_DIPOLE_HINT_RE = re.compile(r"X=.*Y=.*Z=.*Tot=")
_TOT_RE = re.compile(r"^\s*Tot\s+([-\d.DdE+]+)")
_X_RE = re.compile(r"^\s*x\s+([-\d.DdE+]+)")
_Y_RE = re.compile(r"^\s*y\s+([-\d.DdE+]+)")
_Z_RE = re.compile(r"^\s*z\s+([-\d.DdE+]+)")
_ALPHA_STATIC_RE = re.compile(r"^\s*Alpha\(0;0\):")
_ALPHA_DYN_RE = re.compile(r"Alpha\(-w;w\)\s+w=\s*([\d.]+)\s*nm")
_ALPHA_LINE_RE = re.compile(r"^\s*(\w+)\s+([-\d.DdE+]+)")

_UNITS = {
    "au": (1.0, 1.0, "au", "au (bohr^3)"),
    "debye": (AU_TO_DEBYE, 1.0, "Debye", "au (bohr^3)"),
    "si": (8.47836e-30, 1.48185e-25 * 1e6, "10^-30 C*m", "10^-24 esu (cm^3)"),
}


def _parse_alpha_block(lines, k_start, n):
    iso = aniso = xx = yx = yy = zx = zy = zz = math.nan
    for ln in lines[k_start:min(k_start + 10, n)]:
        m = _ALPHA_LINE_RE.match(ln)
        if not m:
            continue
        label = m.group(1).lower()
        val = fortran_to_float(m.group(2))
        if label == "iso":
            iso = val
        elif label == "aniso":
            aniso = val
        elif label == "xx":
            xx = val
        elif label == "yx":
            yx = val
        elif label == "yy":
            yy = val
        elif label == "zx":
            zx = val
        elif label == "zy":
            zy = val
        elif label == "zz":
            zz = val

    if math.isnan(xx):
        return None
    tensor = np.array([[xx, yx, zx], [yx, yy, zy], [zx, zy, zz]])
    return {"iso": iso, "aniso": aniso, "tensor": tensor}


def g16_dipole_polar(filename, units="au", lines=None):
    """Extracts the dipole moment and polarisability from a Gaussian 16
    .out/.log file.

    Parameters
    ----------
    filename : str
    units : {'au' (default), 'debye', 'si'}
    lines : list[str], optional — pre-read file lines (see g16_read_all)

    Returns
    -------
    dp : Struct — mu_x/mu_y/mu_z/mu_tot/mu_units, alpha_iso/alpha_aniso/
        alpha_tensor/alpha_units, alpha_dyn (list of Struct)/N_dyn,
        has_deriv, filename.
    """
    units_key = units.lower()
    if units_key not in _UNITS:
        raise ValueError("g16_dipole_polar: units must be 'au', 'Debye' or 'SI'.")

    if lines is None:
        lines = read_lines(filename)
    n = len(lines)

    mu_xyz_au = [math.nan, math.nan, math.nan]
    mu_tot_au = math.nan
    alpha0_iso = alpha0_aniso = math.nan
    alpha0_tens = np.full((3, 3), math.nan)
    alpha_dyn_list = []

    for i, ln in enumerate(lines):
        if _COMPACT_DIPOLE_HINT_RE.search(ln) and "Dipole moment" in lines[max(0, i - 1)]:
            m = _COMPACT_DIPOLE_RE.search(ln)
            if m:
                mu_xyz_au = [float(m.group(j)) / AU_TO_DEBYE for j in (1, 2, 3)]
                mu_tot_au = float(m.group(4)) / AU_TO_DEBYE
            continue

        hdr1 = lines[max(0, i - 1)]
        hdr2 = lines[max(0, i - 2)]
        is_input_orient = (
            ("Electric dipole moment" in hdr1 and "input orientation" in hdr1)
            or ("Electric dipole moment" in hdr2 and "input orientation" in hdr2)
        )
        if is_input_orient:
            for off in range(5):
                if i + off >= n:
                    break
                ln2 = lines[i + off]
                m = _TOT_RE.match(ln2)
                if m:
                    mu_tot_au = fortran_to_float(m.group(1))
                m = _X_RE.match(ln2)
                if m:
                    mu_xyz_au[0] = fortran_to_float(m.group(1))
                m = _Y_RE.match(ln2)
                if m:
                    mu_xyz_au[1] = fortran_to_float(m.group(1))
                m = _Z_RE.match(ln2)
                if m:
                    mu_xyz_au[2] = fortran_to_float(m.group(1))

        if _ALPHA_STATIC_RE.match(ln):
            block = _parse_alpha_block(lines, i + 1, n)
            if block is not None:
                alpha0_iso = block["iso"]
                alpha0_aniso = block["aniso"]
                alpha0_tens = block["tensor"]
            continue

        m = _ALPHA_DYN_RE.search(ln)
        if m:
            lam_nm = float(m.group(1))
            block = _parse_alpha_block(lines, i + 1, n)
            if block is not None:
                alpha_dyn_list.append(Struct(
                    lambda_nm=lam_nm, freq_au=45.5640 / lam_nm,
                    iso=block["iso"], aniso=block["aniso"], tensor=block["tensor"],
                ))
            continue

    ufac_mu, ufac_alpha, ulbl_mu, ulbl_alpha = _UNITS[units_key]

    for d in alpha_dyn_list:
        d.iso *= ufac_alpha
        d.aniso *= ufac_alpha
        d.tensor = d.tensor * ufac_alpha

    dp = Struct(
        mu_x=mu_xyz_au[0] * ufac_mu, mu_y=mu_xyz_au[1] * ufac_mu,
        mu_z=mu_xyz_au[2] * ufac_mu, mu_tot=mu_tot_au * ufac_mu, mu_units=ulbl_mu,
        alpha_iso=alpha0_iso * ufac_alpha, alpha_aniso=alpha0_aniso * ufac_alpha,
        alpha_tensor=alpha0_tens * ufac_alpha, alpha_units=ulbl_alpha,
        alpha_dyn=alpha_dyn_list, N_dyn=len(alpha_dyn_list),
        has_deriv=False, filename=filename,
    )

    print(f"\n── g16_dipole_polar: {filename} ──")
    print(f"  Dipole  mu = ({dp.mu_x:.4f}, {dp.mu_y:.4f}, {dp.mu_z:.4f})  "
          f"|mu| = {dp.mu_tot:.4f} {dp.mu_units}")
    print(f"  alpha iso   = {dp.alpha_iso:.3f} {dp.alpha_units}")
    print(f"  alpha aniso = {dp.alpha_aniso:.3f} {dp.alpha_units}")
    print("  alpha tensor [au]:")
    for row in alpha0_tens:
        print(f"    {row[0]:10.3f}  {row[1]:10.3f}  {row[2]:10.3f}")
    if dp.N_dyn > 0:
        print("  Alpha(-w;w):")
        for d in dp.alpha_dyn:
            print(f"    lambda={d.lambda_nm:.1f} nm  iso={d.iso:.3f}  aniso={d.aniso:.3f}")
    print()

    return dp
