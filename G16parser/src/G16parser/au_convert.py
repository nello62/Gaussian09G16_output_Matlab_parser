"""Atomic-unit <-> SI (or other common computational-chemistry unit)
conversions. Python port of au2SI/AU_convert.m, au2SI/AU_table.m, and
au2SI/G_gaussian_field_convert.m -- see those files' docstrings for the
full derivation/reference notes. All factors are CODATA 2022 (NIST).

Naming: g16_ prefix on all three, like every other public function in
this package, even though none of the three MATLAB originals have a
G09_/G16_ prefix themselves (they are version-agnostic) -- same
rationale as g16_read_input (see the toolbox manual's own note on that
function): a single Python package needs no version split, so it simply
follows this package's own convention throughout.
"""

import numpy as np

# CODATA 2022 constants (exact where noted)
_e = 1.602176634e-19       # C            (exact)
_me = 9.1093837139e-31     # kg
_mp = 1.67262192595e-27    # kg
_hbar = 1.054571817e-34    # J*s          (exact)
_h = 6.62607015e-34        # J*s          (exact)
_c = 299792458.0           # m/s          (exact)
_kB = 1.380649e-23         # J/K          (exact)
_NA = 6.02214076e23        # mol^-1       (exact)
_a0 = 5.29177210544e-11    # m            (Bohr radius)
_Eh = 4.3597447222060e-18  # J            (Hartree)
_eps0 = 8.8541878188e-12   # F/m


def _energy_factor(target):
    t = target.lower()
    if t in ("kj/mol", "kjmol", "kj"):
        return _Eh * _NA / 1000, "kJ/mol"
    if t in ("kcal/mol", "kcalmol", "kcal"):
        return _Eh * _NA / 4184, "kcal/mol"
    if t in ("ev", "electronvolt"):
        return _Eh / _e, "eV"
    if t in ("cm-1", "cm^-1", "wavenumber", "wn"):
        return _Eh / (_h * _c * 100), "cm-1"
    return _Eh, "J"


def _length_factor(target):
    t = target.lower()
    if t in ("angstrom", "a", "ang", "å"):
        return _a0 * 1e10, "Å"
    if t in ("pm", "picometer"):
        return _a0 * 1e12, "pm"
    return _a0, "m"


def _dipole_factor(target):
    t = target.lower()
    if t in ("debye", "d"):
        return (_e * _a0) / 3.33564095e-30, "Debye"
    return _e * _a0, "C·m"


def _polar_factor(target):
    alpha_si = _e**2 * _a0**2 / _Eh
    t = target.lower()
    if t in ("angstrom3", "a3", "ang3", "å3"):
        return alpha_si / (4 * np.pi * _eps0) * 1e30, "Å³"
    if t in ("esu", "cm3"):
        return alpha_si / (4 * np.pi * _eps0) * 1e6, "cm³ (esu)"
    return alpha_si, "C²·m²·J⁻¹"


def _beta_factor(target):
    beta_si = _e**3 * _a0**3 / _Eh**2
    if target.lower() == "esu":
        return 8.6392e-33, "×10⁻³³ esu"
    return beta_si, "C³·m³·J⁻²"


def _gamma_factor():
    return _e**4 * _a0**4 / _Eh**3, "C⁴·m⁴·J⁻³"


def _force_factor(target):
    F_au = _Eh / _a0
    if target.lower() in ("nn", "nanonewton"):
        return F_au * 1e9, "nN"
    return F_au, "N"


def _frcconst_factor(target):
    k_au = _Eh / _a0**2
    if target.lower() in ("mdyne/a", "mdyne/ang", "mdyne/angstrom", "mdyne"):
        return k_au * 1e-2, "mDyne/Å"
    return k_au, "N/m"


def _frequency_factor(target):
    omega_au = _Eh / _hbar
    t = target.lower()
    if t in ("hz", "hertz"):
        return omega_au / (2 * np.pi), "Hz"
    if t in ("thz", "terahertz"):
        return omega_au / (2 * np.pi) / 1e12, "THz"
    if t in ("cm-1", "cm^-1", "wavenumber", "wn"):
        return omega_au / (2 * np.pi * _c * 100), "cm⁻¹"
    return omega_au, "rad/s"


def _time_factor(target):
    t_au = _hbar / _Eh
    t = target.lower()
    if t in ("fs", "femtosecond"):
        return t_au * 1e15, "fs"
    if t in ("as", "attosecond"):
        return t_au * 1e18, "as"
    return t_au, "s"


def _pressure_factor(target):
    P_au = _Eh / _a0**3
    t = target.lower()
    if t in ("gpa", "gigapascal"):
        return P_au / 1e9, "GPa"
    if t in ("atm", "atmosphere"):
        return P_au / 101325, "atm"
    if t == "bar":
        return P_au / 1e5, "bar"
    return P_au, "Pa"


def _efield_factor(target):
    E_au = _Eh / (_e * _a0)
    t = target.lower()
    if t in ("gv/m", "gvm"):
        return E_au / 1e9, "GV/m"
    if t in ("v/a", "v/angstrom", "v/ang"):
        return E_au * 1e-10, "V/Å"
    if t in ("mv/cm", "mvcm"):
        return E_au / 1e6, "MV/cm"
    return E_au, "V/m"


def _mass_factor(target):
    if target.lower() in ("da", "amu", "u", "dalton"):
        return _me / 1.66053906892e-27, "Da (AMU)"
    return _me, "kg"


_QUANTITY_FACTORS = {
    "energy": lambda t: _energy_factor(t),
    "length": lambda t: _length_factor(t),
    "dipole": lambda t: _dipole_factor(t), "dipolemoment": lambda t: _dipole_factor(t), "mu": lambda t: _dipole_factor(t),
    "polar": lambda t: _polar_factor(t), "polarisability": lambda t: _polar_factor(t), "polarizability": lambda t: _polar_factor(t), "alpha": lambda t: _polar_factor(t),
    "beta": lambda t: _beta_factor(t), "hyperpolar": lambda t: _beta_factor(t), "firsthyperpolar": lambda t: _beta_factor(t),
    "gamma": lambda t: _gamma_factor(), "secondhyperpolar": lambda t: _gamma_factor(),
    "force": lambda t: _force_factor(t),
    "frcconst": lambda t: _frcconst_factor(t), "forceconstant": lambda t: _frcconst_factor(t), "frc": lambda t: _frcconst_factor(t),
    "frequency": lambda t: _frequency_factor(t), "freq": lambda t: _frequency_factor(t), "omega": lambda t: _frequency_factor(t),
    "time": lambda t: _time_factor(t),
    "pressure": lambda t: _pressure_factor(t),
    "efield": lambda t: _efield_factor(t), "electricfield": lambda t: _efield_factor(t),
    "mass": lambda t: _mass_factor(t),
    "charge": lambda t: (_e, "C"),
    "magmom": lambda t: (_e * _hbar / (2 * _me), "J/T"), "magneticmoment": lambda t: (_e * _hbar / (2 * _me), "J/T"),
    "temp": lambda t: (_kB / _Eh, "K (via kBT)"), "temperature": lambda t: (_kB / _Eh, "K (via kBT)"),
}


def g16_au_convert(value, quantity, direction, target="", verbose=None):
    """Converts a physical quantity between atomic units (au) and SI (or
    another common computational-chemistry unit).

    value : scalar or array-like
    quantity : str -- see au_convert.m's docstring for the full list
        ('energy', 'length', 'dipole', 'polar', 'beta', 'gamma', 'force',
        'frcconst', 'frequency', 'time', 'pressure', 'efield', 'mass',
        'charge', 'magmom', 'temp'), or 'help' to print the reference
        table instead of converting anything (returns None).
    direction : 'au2si' (au -> SI/target) | 'si2au' (SI/target -> au)
    target : target unit string (only needed when a quantity has more
        than one SI-side unit; see the quantity's own default above)
    verbose : print a one-line human-readable summary for a scalar
        input (default: True for a scalar value, False for an array,
        matching the MATLAB original's isscalar check)

    Returns the converted value (same shape as value), or None in
    'help' mode.
    """
    if quantity.lower() == "help":
        _print_help()
        return None

    direction = direction.lower().strip()
    if direction not in ("au2si", "si2au"):
        raise ValueError("g16_au_convert: direction must be 'au2si' or 'si2au'.")

    q = quantity.lower().strip()
    if q not in _QUANTITY_FACTORS:
        raise ValueError(
            f"g16_au_convert: unknown quantity '{quantity}'. "
            "Call g16_au_convert(None, 'help', '') for the full list."
        )
    factor, unit_si = _QUANTITY_FACTORS[q](target.lower().strip())

    value_arr = np.asarray(value, dtype=float)
    result = value_arr * factor if direction == "au2si" else value_arr / factor

    is_scalar = np.ndim(value_arr) == 0 or value_arr.size == 1
    if verbose is None:
        verbose = is_scalar
    if verbose:
        v = float(value_arr) if is_scalar else value_arr
        r = float(result) if is_scalar else result
        if direction == "au2si":
            print(f"{v:g} au  ->  {r:g} {unit_si}" if is_scalar else f"[array] au -> {unit_si}")
        else:
            print(f"{v:g} {unit_si}  ->  {r:g} au" if is_scalar else f"[array] {unit_si} -> au")

    return float(result) if (is_scalar and np.ndim(value) == 0) else result


def _print_help():
    print()
    print("=" * 70)
    print("  g16_au_convert -- Atomic Units <-> SI Conversion")
    print("  CODATA 2022 constants  |  Sebastiano Trusso, CNR-IPCF Messina")
    print("=" * 70)
    print()
    print("  Usage: result = g16_au_convert(value, quantity, direction, target=unit)")
    print()
    rows = [
        ("QUANTITY", "au UNIT", "SI/target", "FACTOR (au->SI)"),
        ("-" * 9, "-" * 10, "-" * 20, "-" * 18),
        ("energy", "Hartree", "J", "4.359745e-18"),
        ("", "", "kJ/mol", "2625.4996"),
        ("", "", "kcal/mol", "627.5095"),
        ("", "", "eV", "27.211386"),
        ("", "", "cm-1", "219474.631"),
        ("length", "Bohr (a0)", "m", "5.291772e-11"),
        ("", "", "Angstrom", "0.529177"),
        ("", "", "pm", "52.9177"),
        ("dipole", "e*a0", "C*m", "8.478354e-30"),
        ("", "", "Debye", "2.541747"),
        ("polar", "au", "C^2*m^2*J^-1", "1.648777e-41"),
        ("", "", "Angstrom3", "0.148185"),
        ("", "", "esu (cm3)", "1.48185e-25"),
        ("beta", "au", "C^3*m^3*J^-2", "3.206361e-53"),
        ("", "", "esu", "8.6392e-33"),
        ("gamma", "au", "C^4*m^4*J^-3", "6.235380e-65"),
        ("force", "Eh/a0", "N", "8.238724e-8"),
        ("", "", "nN", "82.38724"),
        ("frcconst", "Eh/a0^2", "N/m", "1556.893"),
        ("", "", "mDyne/A", "15.5689"),
        ("frequency", "Eh/hbar", "rad/s", "4.134137e16"),
        ("", "", "Hz", "6.579684e15"),
        ("", "", "THz", "6579.684"),
        ("", "", "cm-1", "219474.631"),
        ("time", "hbar/Eh", "s", "2.418884e-17"),
        ("", "", "fs", "0.024189"),
        ("", "", "as", "24.1888"),
        ("pressure", "Eh/a0^3", "Pa", "2.942102e13"),
        ("", "", "GPa", "29421.0"),
        ("", "", "atm", "2.904430e8"),
        ("efield", "Eh/(e*a0)", "V/m", "5.142207e11"),
        ("", "", "GV/m", "514.221"),
        ("", "", "V/A", "51.4221"),
        ("mass", "me", "kg", "9.109384e-31"),
        ("", "", "Da (AMU)", "5.485799e-4"),
        ("charge", "e", "C", "1.602177e-19"),
        ("magmom", "muB", "J/T", "9.274010e-24"),
        ("temp", "Eh", "K  (via kBT)", "3.157750e5"),
    ]
    for r in rows:
        print(f"  {r[0]:<14} {r[1]:<12} {r[2]:<22} {r[3]}")
    print()
    print("Examples:")
    print("  g16_au_convert(-875.932, 'energy',    'au2si', target='kJ/mol')")
    print("  g16_au_convert(6.3347,   'dipole',    'si2au', target='Debye')")
    print("  g16_au_convert(259.03,   'polar',     'au2si', target='Angstrom3')")
    print("  g16_au_convert(1582.8,   'frequency', 'si2au', target='cm-1')")
    print("  g16_au_convert(mol.xyz,  'length',    'au2si', target='Angstrom')")
    print("  g16_au_convert(None,     'help',      '')")
    print()


def g16_au_table(section=""):
    """Prints a formatted reference table of atomic-unit conversion
    factors (CODATA 2022).

    section : '' (default, prints every section) | one of the
        g16_au_convert quantity strings (e.g. 'energy', 'length') to
        print only that section | 'constants' for the raw fundamental
        constants | 'all' (same as '').
    """
    filt = section.lower().strip()

    alpha_au = _e**2 * _a0**2 / _Eh
    beta_au = _e**3 * _a0**3 / _Eh**2
    gamma_au = _e**4 * _a0**4 / _Eh**3

    sections = [
        ("FUNDAMENTAL CONSTANTS  (CODATA 2022)", "constants", [
            ("Bohr radius  a0", f"{_a0:.10e} m"),
            ("Hartree  Eh", f"{_Eh:.10e} J"),
            ("Electron mass  me", f"{_me:.10e} kg"),
            ("Proton mass  mp", f"{_mp:.10e} kg"),
            ("Elementary charge  e", f"{_e:.10e} C"),
            ("Planck const  h", f"{_h:.10e} J*s"),
            ("Reduced Planck  hbar", f"{_hbar:.10e} J*s"),
            ("Speed of light  c", f"{int(_c)} m/s (exact)"),
            ("Boltzmann  kB", f"{_kB:.10e} J/K (exact)"),
            ("Avogadro  NA", f"{_NA:.10e} mol^-1 (exact)"),
            ("Vacuum permittivity eps0", f"{_eps0:.10e} F/m"),
        ]),
        ("ENERGY", "energy", [
            ("1 Hartree (Eh)", f"= {_Eh:.10e} J"),
            ("", f"= {_Eh*_NA/1000:.6f} kJ/mol"),
            ("", f"= {_Eh*_NA/4184:.6f} kcal/mol"),
            ("", f"= {_Eh/_e:.6f} eV"),
            ("", f"= {_Eh/(_h*_c*100):.4f} cm^-1"),
            ("", f"= {_Eh/_e*1000:.6f} meV"),
            ("1 eV", f"= {_e:.10e} J"),
            ("", f"= {_e*_NA/1000:.4f} kJ/mol"),
            ("", f"= {_e/(_h*_c*100):.4f} cm^-1"),
            ("1 cm^-1", f"= {_h*_c*100:.10e} J"),
            ("", f"= {_h*_c*100/_Eh:.8e} Hartree"),
            ("1 kJ/mol", f"= {1000/(_Eh*_NA):.8e} Hartree"),
            ("1 kcal/mol", f"= {4184/(_Eh*_NA):.8e} Hartree"),
            ("1 kBT (298.15K)", f"= {_kB*298.15/_Eh:.8e} Hartree  (= {_kB*298.15*_NA/1000:.4f} kJ/mol)"),
        ]),
        ("LENGTH", "length", [
            ("1 Bohr (a0)", f"= {_a0:.10e} m"),
            ("", f"= {_a0*1e10:.10f} A"),
            ("", f"= {_a0*1e12:.6f} pm"),
            ("1 A", f"= {1e-10/_a0:.10f} Bohr"),
            ("1 pm", f"= {1e-12/_a0:.10f} Bohr"),
            ("1 nm", f"= {1e-9/_a0:.6f} Bohr"),
        ]),
        ("DIPOLE MOMENT", "dipole", [
            ("1 au  (e*a0)", f"= {_e*_a0:.10e} C*m"),
            ("", f"= {_e*_a0/3.33564095e-30:.8f} Debye"),
            ("1 Debye", f"= {3.33564095e-30:.10e} C*m"),
            ("", f"= {3.33564095e-30/(_e*_a0):.8f} au"),
        ]),
        ("POLARISABILITY  alpha", "polar", [
            ("1 au  (e^2 a0^2/Eh)", f"= {alpha_au:.10e} C^2*m^2*J^-1"),
            ("", f"= {alpha_au/(4*np.pi*_eps0)*1e30:.6f} A^3  (volume)"),
            ("", f"= {alpha_au/(4*np.pi*_eps0)*1e6:.6e} cm^3  (esu)"),
            ("1 A^3", f"= {(4*np.pi*_eps0)*1e-30/alpha_au:.6f} au"),
        ]),
        ("FIRST HYPERPOLARISABILITY  beta", "beta", [
            ("1 au  (e^3 a0^3/Eh^2)", f"= {beta_au:.6e} C^3*m^3*J^-2  (SI)"),
            ("", "= 8.6392e-33 esu  (x1e-33)"),
            ("1 x1e-33 esu", f"= {1e-33/8.6392e-33:.6f} au"),
        ]),
        ("SECOND HYPERPOLARISABILITY  gamma", "gamma", [
            ("1 au  (e^4 a0^4/Eh^3)", f"= {gamma_au:.6e} C^4*m^4*J^-3  (SI)"),
        ]),
        ("FORCE", "force", [
            ("1 au  (Eh/a0)", f"= {_Eh/_a0:.6e} N"),
            ("", f"= {_Eh/_a0*1e9:.6f} nN"),
            ("1 nN", f"= {1e-9/(_Eh/_a0):.8f} au"),
        ]),
        ("FORCE CONSTANT", "frcconst", [
            ("1 au  (Eh/a0^2)", f"= {_Eh/_a0**2:.4f} N/m"),
            ("", f"= {_Eh/_a0**2*1e-2:.4f} mDyne/A"),
            ("1 mDyne/A", f"= {1e-2/(_Eh/_a0**2):.8f} au"),
            ("1 N/m", f"= {1/(_Eh/_a0**2):.8f} au"),
        ]),
        ("FREQUENCY", "frequency", [
            ("1 au  (Eh/hbar)", f"= {_Eh/_hbar:.6e} rad/s"),
            ("", f"= {_Eh/(_hbar*2*np.pi):.6e} Hz"),
            ("", f"= {_Eh/(_hbar*2*np.pi*_c*100):.4f} cm^-1"),
            ("", f"= {_Eh/(_hbar*2*np.pi*1e12):.4f} THz"),
            ("1 cm^-1", f"= {_h*_c*100/_Eh:.10e} au"),
            ("1 THz", f"= {1e12*2*np.pi*_hbar/_Eh:.10e} au"),
        ]),
        ("TIME", "time", [
            ("1 au  (hbar/Eh)", f"= {_hbar/_Eh:.6e} s"),
            ("", f"= {_hbar/_Eh*1e18:.6f} as  (attoseconds)"),
            ("", f"= {_hbar/_Eh*1e15:.6e} fs"),
            ("1 fs", f"= {1e-15/(_hbar/_Eh):.4f} au"),
            ("1 as", f"= {1e-18/(_hbar/_Eh):.6f} au"),
        ]),
        ("PRESSURE", "pressure", [
            ("1 au  (Eh/a0^3)", f"= {_Eh/_a0**3:.6e} Pa"),
            ("", f"= {_Eh/_a0**3/1e9:.4f} GPa"),
            ("", f"= {_Eh/_a0**3/101325:.4e} atm"),
            ("1 GPa", f"= {1e9/(_Eh/_a0**3):.6e} au"),
        ]),
        ("ELECTRIC FIELD", "efield", [
            ("1 au  (Eh/ea0)", f"= {_Eh/(_e*_a0):.6e} V/m"),
            ("", f"= {_Eh/(_e*_a0)/1e9:.4f} GV/m"),
            ("", f"= {_Eh/(_e*_a0)*1e-10:.4f} V/A"),
            ("1 V/A", f"= {1e10*_e*_a0/_Eh:.6e} au"),
        ]),
        ("MASS", "mass", [
            ("1 au  (me)", f"= {_me:.10e} kg"),
            ("", f"= {_me/1.66053906892e-27:.8e} Da (AMU)"),
            ("1 Da (AMU)", f"= {1.66053906892e-27/_me:.4f} au  (me)"),
            ("1 mp", f"= {_mp/_me:.4f} au  (me)"),
        ]),
        ("CHARGE", "charge", [
            ("1 au  (e)", f"= {_e:.10e} C  (elementary charge)"),
        ]),
        ("MAGNETIC MOMENT", "magmom", [
            ("1 au  (muB)", f"= {_e*_hbar/(2*_me):.10e} J/T  (Bohr magneton)"),
            ("1 muN", f"= {_e*_hbar/(2*_mp):.6e} J/T  (nuclear magneton)"),
        ]),
        ("TEMPERATURE  (via thermal energy)", "temp", [
            ("1 Eh = kBT  <->  T", f"= {_Eh/_kB:.4f} K"),
            ("1 K  <->  kBT", f"= {_kB/_Eh:.6e} Eh"),
            ("kBT at 298.15 K", f"= {_kB*298.15/_Eh:.6e} Eh  = {_kB*298.15/_e*1000:.4f} meV"),
        ]),
    ]

    print()
    print("=" * 70)
    print("  g16_au_table -- Atomic Unit Conversion Reference  (CODATA 2022)")
    print("=" * 70)
    for title, tag, rows in sections:
        if filt and filt != tag and filt != "all":
            continue
        print(f"\n  -- {title} --")
        for lbl, val in rows:
            if lbl:
                print(f"  {lbl:<28}  {val}")
            else:
                print(f"  {'':<28}  {val}")
    print("\n  For conversions: g16_au_convert(value, quantity, 'au2si'|'si2au', target=unit)")
    print("  For help:        g16_au_convert(None, 'help', '')\n")


def g16_gaussian_field_convert(value, direction, unit="au", verbose=True):
    """Converts between Gaussian's Field=X+N route keyword integer and a
    physical electric-field strength (V/Angstrom, V/m, V/cm, or a.u.).

    Gaussian's Field keyword (e.g. Field=X+10) specifies a static field
    whose magnitude, in atomic units, is N*0.0001 (verified against
    gaussian.com/field/).

    direction : 'g2phys' -- value is the Gaussian keyword integer N;
        returns the field strength in `unit` (default 'au').
                'phys2g' -- value is a field strength in `unit`; returns
        the Gaussian keyword integer N, rounded to the nearest integer
        (Gaussian's keyword can only express multiples of 0.0001 a.u.).
        A warning is printed if rounding changes the value by > 0.5%.
    unit : 'au' | 'V/m' | 'V/cm' | 'V/Ang' (default 'au')
    verbose : print a one-line summary (default True)

    value may be a numpy array; the conversion is applied elementwise.

    Example
    -------
    g16_gaussian_field_convert(10, 'g2phys')                     # 0.0010 (au)
    g16_gaussian_field_convert(10, 'g2phys', unit='V/Ang')       # 0.5142 V/Ang
    g16_gaussian_field_convert(0.0025, 'phys2g', unit='au')      # 25
    """
    direction = direction.lower()
    if direction not in ("g2phys", "phys2g"):
        raise ValueError("g16_gaussian_field_convert: direction must be 'g2phys' or 'phys2g'.")
    if unit.lower() not in ("au", "v/m", "v/cm", "v/ang"):
        raise ValueError("g16_gaussian_field_convert: unit must be 'au', 'V/m', 'V/cm', or 'V/Ang'.")

    au2Vm = _Eh / (_e * _a0)
    u = unit.lower()
    if u == "au":
        unit_per_au, unit_label = 1.0, "au"
    elif u == "v/m":
        unit_per_au, unit_label = au2Vm, "V/m"
    elif u == "v/cm":
        unit_per_au, unit_label = au2Vm / 100, "V/cm"
    else:  # v/ang
        unit_per_au, unit_label = au2Vm * 1e-10, "V/Ang"

    value_arr = np.asarray(value, dtype=float)

    if direction == "g2phys":
        field_au = value_arr * 1e-4
        out = field_au * unit_per_au
        if verbose:
            print(f"g16_gaussian_field_convert: Field keyword N={value_arr}  ->  {out} {unit_label}")
    else:
        field_au = value_arr / unit_per_au
        N_raw = field_au / 1e-4
        out = np.round(N_raw)
        rel_err = np.abs(out - N_raw) / np.maximum(np.abs(N_raw), np.finfo(float).eps)
        if verbose:
            print(f"g16_gaussian_field_convert: {value_arr} {unit_label}  ->  Field keyword N={out}")
        if np.any(rel_err > 0.005):
            import warnings
            warnings.warn(
                "The requested field strength is not an exact multiple of 0.0001 a.u. "
                "(Gaussian's Field keyword granularity); rounded to the nearest integer N, "
                f"a relative change of up to {100*np.max(rel_err):.2f}%."
            )

    is_scalar = np.ndim(value_arr) == 0
    return (int(out) if direction == "phys2g" else float(out)) if is_scalar else out
