import math
import re

from ._common import Struct, read_lines, fortran_to_float

_VIB_HDR = "Diagonal vibrational hyperpolarizability:"
_BETA0_RE = re.compile(r"^\s*Beta\(0;0,0\):")
_BETA_DYN_RE = re.compile(r"Beta\(-w;w,0\)\s+w=\s*([\d.]+)\s*nm")
_LABEL_RE = re.compile(r"^\s*(\S+)\s+([-\d.DdE+]+)")

_TENSOR_FIELDS = ["xxx", "xxy", "yxy", "yyy", "xxz", "yxz", "yyz", "zxz", "zyz", "zzz",
                  "yxx", "zxx", "zyx", "zzx", "zyy", "zzy"]

_UNITS = {
    "au": (1.0, "au"),
    "esu": (8.6392e-3, "10^-30 esu"),
    "si": (3.2063e-3, "10^-50 C^3m^3J^-2"),
}


def _parse_beta_block(lines, k_start, n):
    par_z = perp_z = vec_x = vec_y = vec_z = math.nan
    tensor = {f: math.nan for f in _TENSOR_FIELDS}

    for ln in lines[k_start:min(k_start + 25, n)]:
        if ln.strip() == "":
            break
        if re.match(r"^\s*-{10,}", ln):
            break
        m = _LABEL_RE.match(ln)
        if not m:
            continue
        label, val_s = m.group(1), m.group(2)
        val = fortran_to_float(val_s)
        if label == "||(z)":
            par_z = val
        elif label == "_|_(z)":
            perp_z = val
        elif label == "x":
            vec_x = val
        elif label == "y":
            vec_y = val
        elif label == "z":
            vec_z = val
        elif label in tensor:
            tensor[label] = val

    if math.isnan(vec_x):
        return None

    beta_vec = math.sqrt(vec_x ** 2 + vec_y ** 2 + vec_z ** 2)
    return {
        "par_z": par_z, "perp_z": perp_z, "vec_x": vec_x, "vec_y": vec_y, "vec_z": vec_z,
        "beta_vec": beta_vec, "tensor": tensor, "units": "au",
    }


def _apply_units(raw, ufac, ulbl):
    if raw is None:
        return Struct(par_z=math.nan, perp_z=math.nan, vec_x=math.nan, vec_y=math.nan,
                       vec_z=math.nan, beta_vec=math.nan, tensor=None, units=ulbl)
    tensor = {k: (v * ufac if not math.isnan(v) else v) for k, v in raw["tensor"].items()}
    return Struct(
        par_z=raw["par_z"] * ufac, perp_z=raw["perp_z"] * ufac,
        vec_x=raw["vec_x"] * ufac, vec_y=raw["vec_y"] * ufac, vec_z=raw["vec_z"] * ufac,
        beta_vec=raw["beta_vec"] * ufac, tensor=Struct(**tensor), units=ulbl,
    )


def g16_hyperpolar(filename, units="au", lines=None):
    """Extracts the dipole hyperpolarisability Beta from a Gaussian 16
    .out/.log file (static Beta(0;0,0), dynamic Beta(-w;w,0) per laser
    wavelength, and the vibrational contribution if present).

    Returns
    -------
    hp : Struct — beta0 (Struct), beta_dyn (list[Struct]), N_dyn,
        beta_vib (list[float] or None), has_vib, filename.
    """
    units_key = units.lower()
    if units_key not in _UNITS:
        raise ValueError("g16_hyperpolar: units must be 'au', 'esu' or 'SI'.")
    ufac, ulbl = _UNITS[units_key]

    if lines is None:
        lines = read_lines(filename)
    n = len(lines)

    beta0_raw = None
    beta_dyn_raw = []
    beta_vib = None

    i = 0
    while i < n:
        ln = lines[i]

        if _VIB_HDR in ln:
            if i + 1 < n:
                vals = [float(x) for x in lines[i + 1].split()]
                if len(vals) >= 3:
                    beta_vib = vals[:3]
            i += 2
            continue

        if _BETA0_RE.match(ln):
            beta0_raw = _parse_beta_block(lines, i + 1, n)
            i += 15
            continue

        m = _BETA_DYN_RE.search(ln)
        if m:
            lam = float(m.group(1))
            bdata = _parse_beta_block(lines, i + 1, n)
            if bdata is not None:
                beta_dyn_raw.append((lam, bdata))
            i += 15
            continue

        i += 1

    hp = Struct(
        beta0=_apply_units(beta0_raw, ufac, ulbl),
        beta_dyn=[], N_dyn=len(beta_dyn_raw),
        beta_vib=None, has_vib=False, filename=filename,
    )
    for lam, bdata in beta_dyn_raw:
        entry = _apply_units(bdata, ufac, ulbl)
        entry.lambda_nm = lam
        hp.beta_dyn.append(entry)

    if beta_vib is not None:
        hp.beta_vib = [v * ufac for v in beta_vib]
        hp.has_vib = True

    print(f"\n── g16_hyperpolar: {filename} ──")
    if beta0_raw is not None:
        print(f"  Beta(0;0,0):  |beta_vec| = {hp.beta0.beta_vec:.2f} {ulbl}")
        print(f"    || (z) = {hp.beta0.par_z:.2f}   _|_(z) = {hp.beta0.perp_z:.2f}")
    for entry in hp.beta_dyn:
        print(f"  Beta(-w;w,0) lambda={entry.lambda_nm:.1f} nm:  "
              f"|beta_vec| = {entry.beta_vec:.2f} {ulbl}")
    if hp.has_vib:
        print(f"  beta_vib (diag) = [{hp.beta_vib[0]:.2f}  {hp.beta_vib[1]:.2f}  "
              f"{hp.beta_vib[2]:.2f}] {ulbl}")
    print()

    return hp
