from ._common import Struct, read_lines
from .charges import g16_charges
from .energy import g16_energy
from .structure import g16_structure
from .dipole_polar import g16_dipole_polar
from .nmodes import g16_nmodes
from .spectra import g16_spectra
from .route import g16_route
from .charge_mult import g16_charge_mult


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
        .nmodes     — vibrational normal modes (g16_nmodes)
        .spectra    — IR/Raman spectra (g16_spectra)
        .route      — Gaussian route section (g16_route)
        .chargemol  — Struct(charge=, mol=) from g16_charge_mult
    """
    lines = read_lines(filename)

    charge = g16_charges(filename, plot=False, lines=lines)
    energy = g16_energy(filename, lines=lines)
    structure = g16_structure(filename, lines=lines)
    dipolar = g16_dipole_polar(filename, lines=lines)
    nmodes = g16_nmodes(filename, lines=lines)
    spectra = g16_spectra(filename, lines=lines)
    route = g16_route(filename, lines=lines)
    c, m = g16_charge_mult(filename, lines=lines)

    return Struct(
        charge=charge, energy=energy, structure=structure, dipolar=dipolar,
        nmodes=nmodes, spectra=spectra, route=route,
        chargemol=Struct(charge=c, mol=m),
    )
