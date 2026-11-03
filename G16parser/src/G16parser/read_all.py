import warnings

from ._common import Struct, read_lines
from .charges import g16_charges
from .energy import g16_energy
from .structure import g16_structure
from .dipole_polar import g16_dipole_polar
from .nmodes import g16_nmodes
from .spectra import g16_spectra
from .route import g16_route
from .charge_mult import g16_charge_mult
from ._pathjobs import split_jobs


def g16_read_all(filename):
    """Collects all relevant data from a Gaussian 16 .out file, running the
    full set of g16_* extraction functions and assembling the results into
    a single Struct — a one-call summary of a calculation.

    The file is read from disk once and the parsed lines are reused by
    every sub-function (via their `lines=` parameter), instead of each one
    re-reading and re-splitting the file independently.

    Returns
    -------
    T : Struct with fields:
        .charge     — Mulliken and APT atomic charges (g16_charges; plotting disabled)
        .energy     — SCF energy and thermochemistry (g16_energy)
        .structure  — optimised molecular structure (g16_structure)
        .dipolar    — dipole moment and polarisability (g16_dipole_polar)
        .nmodes     — vibrational normal modes (g16_nmodes); None if the
            file has no frequency calculation (e.g. a TD-DFT-only single
            point job with no 'freq' keyword) — a non-blocking warning
            is issued in that case
        .spectra    — IR/Raman spectra (g16_spectra); None under the
            same condition as .nmodes above
        .route      — Gaussian route section (g16_route)
        .chargemol  — Struct(charge=, mol=) from g16_charge_mult

    Notes
    -----
    Every field is taken from the LAST occurrence in the file. For files
    holding several structures -- a relaxed scan, an IRC, or several jobs
    chained with --Link1-- -- the fields may come from different jobs or
    points; a warning is issued, pointing to g16_scan_read / g16_irc_read.
    """
    lines = read_lines(filename)
    _warn_multistructure(filename, lines)

    charge = g16_charges(filename, plot=False, lines=lines)
    energy = g16_energy(filename, lines=lines)
    structure = g16_structure(filename, lines=lines)
    dipolar = g16_dipole_polar(filename, lines=lines)

    # Vibrational normal modes, if present (requires 'freq' in the route
    # section, e.g. absent for a TD-DFT-only single point job) -- caught
    # here rather than left to crash the whole call, since write_report
    # already expects nmodes to be optional (getattr(...) is not None check)
    try:
        nmodes = g16_nmodes(filename, lines=lines)
    except ValueError as exc:
        warnings.warn(f"No vibrational normal modes found in {filename} ({exc}); nmodes omitted.")
        nmodes = None

    # IR/Raman spectra -- same optional-section handling as nmodes above
    try:
        spectra = g16_spectra(filename, lines=lines)
    except ValueError as exc:
        warnings.warn(f"No IR/Raman spectra found in {filename} ({exc}); spectra omitted.")
        spectra = None

    route = g16_route(filename, lines=lines)
    c, m = g16_charge_mult(filename, lines=lines)

    return Struct(
        charge=charge, energy=energy, structure=structure, dipolar=dipolar,
        nmodes=nmodes, spectra=spectra, route=route,
        chargemol=Struct(charge=c, mol=m),
    )


def _warn_multistructure(filename, lines):
    """Warns when g16_read_all's 'last occurrence' logic would mix data
    from several jobs or path points."""
    n_jobs = len(split_jobs(lines))
    has_scan = any("on scan point" in ln for ln in lines)
    has_irc = any("IRC-IRC-IRC" in ln for ln in lines)
    if has_scan:
        warnings.warn(f"{filename} contains a relaxed PES scan: g16_read_all returns only the "
                      "last scan point. Use g16_scan_read for the full scan.")
    if has_irc:
        warnings.warn(f"{filename} contains an IRC: g16_read_all fields may mix the TS "
                      "(frequencies, thermochemistry) with the last IRC point (energy, "
                      "structure). Use g16_irc_read for the reaction path.")
    if n_jobs > 1:
        warnings.warn(f"{filename} contains {n_jobs} jobs chained with --Link1--: g16_read_all "
                      "takes each field from its last occurrence, which may belong to "
                      "different jobs. Split the file (one job per file) or use "
                      "g16_scan_read / g16_irc_read.")
