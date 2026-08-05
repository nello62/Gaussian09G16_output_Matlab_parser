import pandas as pd

from .route import g16_route

_MARK = {True: "YES", False: "no "}


def _has_kw(tokens, name):
    return any(t == name or t.startswith(name + "=") or t.startswith(name + "(") for t in tokens)


def g16_available_data(filename):
    """Lists which toolbox quantities/functions are expected to be
    extractable from a Gaussian 16 output file, based on the route
    section keywords actually present.

    Reads the route section (see g16_route) and checks it against the
    keyword requirements documented in the toolbox manual ("Which .in
    keywords generate which .out section"), so you can tell up front
    which data will actually be present in the file -- e.g. a TD-DFT-only
    single point job (td=(singlets,nstates=100), no 'freq') has no
    vibrational normal modes, no IR/Raman spectra, and no thermochemistry
    (ZPE/H/G/S), so calling g16_nmodes/g16_spectra on it raises a
    ValueError and g16_energy's *_corr/E0/U/H/G/T/P/S fields come back
    NaN -- none of that is a bug, it is simply not in the file. This
    function lets you check that in one call before running the
    extraction functions, rather than discovering it one error/NaN at a
    time.

    This predicts availability from the route keywords alone; it does
    not itself read the data sections, so an incomplete/crashed job can
    still show a keyword as present without the corresponding output
    actually having been printed (the extraction function will still
    raise its own clear error in that case).

    Returns
    -------
    T : pandas.DataFrame with columns Function, Available (bool),
        Requires, Notes. Also prints a formatted list to stdout.

    Example
    -------
        g16.g16_available_data('V_E00_R_24_TD.out')
        T = g16.g16_available_data('molecule.out')
        T[~T['Available']]   # see what's NOT available in this file
    """
    route = g16_route(filename)
    tokens = route.lower().split()

    has_freq = _has_kw(tokens, "freq")
    has_opt = _has_kw(tokens, "opt")
    has_polar = _has_kw(tokens, "polar")
    has_td = _has_kw(tokens, "td")
    has_freq_raman = any(
        (t.startswith("freq=") or t.startswith("freq(")) and "raman" in t for t in tokens
    )
    has_cphf_rdfreq = any(t.startswith("cphf=") and "rdfreq" in t for t in tokens)
    has_field = any(t.startswith("field=") for t in tokens)
    has_nbo = any(t.startswith("pop=") and "nbo" in t for t in tokens)

    rows = [
        ("structure", True, "(default)", "nosymm/nosym suppresses Standard orientation only"),
        ("energy (SCF)", True, "(default)", ""),
        ("energy (thermochemistry: ZPE/H/G/S)", has_freq, "freq", ""),
        ("charges (Mulliken)", True, "(default)", ""),
        ("charges (APT)", has_freq, "freq", "printed automatically with any freq job"),
        ("dipole_polar (dipole moment)", True, "(default)", ""),
        ("dipole_polar (static Alpha)", has_polar, "polar", ""),
        ("dipole_polar (dynamic Alpha)", has_polar and has_cphf_rdfreq, "polar + cphf=rdfreq",
         "also needs an explicit frequency list in the .in file"),
        ("nmodes", has_freq, "freq", ""),
        ("spectra (IR)", has_freq, "freq", ""),
        ("spectra (Raman)", has_freq_raman, "freq=raman", ""),
        ("orbital_energies", True, "(default)", "full eigenvalue list printed unconditionally"),
        ("convergence", has_opt, "opt", ""),
        ("charge_mult", True, "(default)", ""),
        ("route", True, "(default)", ""),
        ("get_bond_length", True, "(default, derived)", "purely geometric, not a real Gaussian section"),
        ("nbo_bonds", has_nbo, "pop=nbo", ""),
        ("gaussian_version", True, "(default)", ""),
        ("restart", True, "(any completed step)", ""),
        ("hyperpolar (vibrational Beta)", has_polar and has_freq, "polar + freq", ""),
        ("hyperpolar (electronic Beta)", has_polar and has_field and has_cphf_rdfreq,
         "polar + field=... + cphf=rdfreq", ""),
        ("tddft", has_td, "td", ""),
    ]

    T = pd.DataFrame(rows, columns=["Function", "Available", "Requires", "Notes"])

    print(f"\n-- g16_available_data: {filename} --")
    print(f"  Route: {route}\n")
    name_w = T["Function"].str.len().max()
    req_w = T["Requires"].str.len().max()
    for _, row in T.iterrows():
        print(f"  {_MARK[row['Available']]}  {row['Function']:<{name_w}}  "
              f"needs: {row['Requires']:<{req_w}}  {row['Notes']}")
    print()

    return T
