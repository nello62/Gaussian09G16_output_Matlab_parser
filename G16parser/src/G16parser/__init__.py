"""G16 Python toolbox — parses and visualises Gaussian 16 .out/.log files.

Python port of the G16/ MATLAB toolbox: data extraction, static
matplotlib plots, and an interactive Tkinter vibrational-mode viewer.

    from G16parser import g16_structure, g16_charges, g16_read_all
    mol = g16_structure('molecule.out')
    T = g16_read_all('molecule.out')
"""

from ._common import Struct
from .gaussian_version import g16_gaussian_version
from .charge_mult import g16_charge_mult
from .route import g16_route
from .structure import g16_structure
from .read_input import g16_read_input
from .energy import g16_energy
from .convergence import g16_convergence
from .dipole_polar import g16_dipole_polar
from .charges import g16_charges
from .nmodes import g16_nmodes
from .spectra import g16_spectra
from .orbital_energies import g16_orbital_energies
from .get_bond_length import g16_get_bond_length
from .hyperpolar import g16_hyperpolar
from .tddft import g16_tddft
from .read_all import g16_read_all
from .draw_molecule import g16_draw_molecule
from .draw_mode import g16_draw_mode
from .draw_orbital import g16_draw_orbital
from .animate_mode import g16_animate_mode
from .mode_viewer import g16_mode_viewer
from .toolbox_list import g16_list
from .write_report import g16_write_report

__all__ = [
    "Struct",
    "g16_gaussian_version", "g16_charge_mult", "g16_route",
    "g16_structure", "g16_read_input", "g16_energy", "g16_convergence", "g16_dipole_polar",
    "g16_charges", "g16_nmodes", "g16_spectra", "g16_orbital_energies",
    "g16_get_bond_length", "g16_hyperpolar", "g16_tddft", "g16_read_all",
    "g16_draw_molecule", "g16_draw_mode", "g16_draw_orbital",
    "g16_animate_mode", "g16_mode_viewer", "g16_list", "g16_write_report",
]
