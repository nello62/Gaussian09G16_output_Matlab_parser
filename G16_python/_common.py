"""Shared helpers used across the G16 Python toolbox.

Not part of the public API — import functions from the top-level
``G16_python`` package instead (e.g. ``from G16_python import g16_structure``).
"""

import re

# Registers the '3d' projection with matplotlib. Required explicitly on
# older matplotlib releases (<3.2ish) where it isn't auto-registered just
# by passing projection="3d" to add_subplot.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


class Struct:
    """A MATLAB-struct-like container: attribute access (``obj.xyz``) backed
    by a plain dict, so results can also be inspected/iterated like one
    (``vars(obj)``, ``obj.__dict__``).
    """

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __repr__(self):
        fields = ", ".join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"Struct({fields})"

    def __eq__(self, other):
        return isinstance(other, Struct) and self.__dict__ == other.__dict__


def read_lines(filename):
    """Reads a Gaussian 16 .out/.log file and returns a list of lines
    (line-ending stripped), matching G16_*.m's fread+strsplit convention.
    """
    with open(filename, "r", encoding="utf-8", errors="replace") as f:
        raw = f.read()
    raw = raw.replace("\r\n", "\n").replace("\r", "\n")
    return raw.split("\n")


def fortran_to_float(s):
    """Converts Fortran D-exponent notation ('0.123D+02') to a Python float."""
    return float(s.replace("D", "E").replace("d", "e"))


# ---------------------------------------------------------------------------
# Atomic number -> symbol table (Z = 1..118)
# ---------------------------------------------------------------------------
_SYMBOLS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]


def z_to_symbol(z):
    return _SYMBOLS[int(z) - 1]


# ---------------------------------------------------------------------------
# CPK colours and covalent radii (Angstrom), shared by draw_molecule and
# get_bond_length — same tables as G16_draw_molecule.m / G16_get_bond_length.m
# ---------------------------------------------------------------------------
CPK_COLORS = {
    "H": (0.60, 0.80, 1.00), "C": (0.30, 0.30, 0.30), "N": (0.10, 0.30, 0.90),
    "O": (0.90, 0.10, 0.10), "F": (0.20, 0.80, 0.20), "P": (1.00, 0.50, 0.00),
    "S": (1.00, 0.85, 0.00), "Cl": (0.20, 0.85, 0.20), "Br": (0.55, 0.20, 0.10),
    "I": (0.45, 0.00, 0.65), "Au": (1.00, 0.82, 0.14), "Ag": (0.75, 0.75, 0.75),
    "Fe": (0.80, 0.40, 0.00), "Zn": (0.50, 0.70, 0.50), "Ca": (0.60, 0.60, 0.60),
    "Mg": (0.50, 0.80, 0.20), "Na": (0.65, 0.40, 0.90), "K": (0.55, 0.20, 0.85),
    "Si": (0.60, 0.50, 0.40), "B": (1.00, 0.65, 0.50), "Cu": (0.72, 0.45, 0.20),
}
DEFAULT_COLOR = (0.65, 0.20, 0.80)

COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "S": 1.05,
    "P": 1.07, "Cl": 1.02, "Br": 1.20, "I": 1.39, "Si": 1.11, "B": 0.84,
    "Na": 1.66, "K": 2.03, "Mg": 1.41, "Ca": 1.76, "Au": 1.36, "Ag": 1.45,
    "Cu": 1.32, "Zn": 1.22, "Cd": 1.44, "Pt": 1.36, "Pd": 1.39, "Ni": 1.24,
    "Fe": 1.32, "Co": 1.26,
}
DEFAULT_RADIUS = 0.80


def get_color(symbol):
    return CPK_COLORS.get(symbol, DEFAULT_COLOR)


def get_radius(symbol):
    return COVALENT_RADII.get(symbol, DEFAULT_RADIUS)


# ---------------------------------------------------------------------------
# Unit conversions
# ---------------------------------------------------------------------------
HARTREE_TO_EV = 27.211386245988
HARTREE_TO_KJMOL = 2625.4996
AU_TO_DEBYE = 2.541746
DEBYE_TO_AU = 1.0 / AU_TO_DEBYE


_FLOAT_RE = re.compile(r"-?\d+\.\d+")


def extract_floats(line):
    """All signed decimal numbers in a line, e.g. eigenvalue rows."""
    return [float(x) for x in _FLOAT_RE.findall(line)]
