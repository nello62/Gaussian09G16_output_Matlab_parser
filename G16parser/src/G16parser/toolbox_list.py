from .animate_mode import g16_animate_mode
from .charge_mult import g16_charge_mult
from .charges import g16_charges
from .convergence import g16_convergence
from .dipole_polar import g16_dipole_polar
from .draw_mode import g16_draw_mode
from .draw_molecule import g16_draw_molecule
from .draw_orbital import g16_draw_orbital
from .energy import g16_energy
from .gaussian_version import g16_gaussian_version
from .get_bond_length import g16_get_bond_length
from .hyperpolar import g16_hyperpolar
from .mode_viewer import g16_mode_viewer
from .nmodes import g16_nmodes
from .orbital_energies import g16_orbital_energies
from .read_all import g16_read_all
from .route import g16_route
from .spectra import g16_spectra
from .structure import g16_structure
from .tddft import g16_tddft
from .write_report import g16_write_report

_FUNCTIONS = [
    g16_gaussian_version, g16_charge_mult, g16_route, g16_structure, g16_energy,
    g16_convergence, g16_dipole_polar, g16_charges, g16_nmodes, g16_spectra,
    g16_orbital_energies, g16_get_bond_length, g16_hyperpolar, g16_tddft,
    g16_read_all, g16_draw_molecule, g16_draw_mode, g16_draw_orbital,
    g16_animate_mode, g16_mode_viewer, g16_write_report,
]


def _description(fn):
    doc = fn.__doc__
    if not doc:
        return ""
    first_para = doc.strip().split("\n\n", 1)[0]
    return " ".join(line.strip() for line in first_para.splitlines())


def g16_list():
    """Lists every function in the G16parser toolbox with a one-line description.

    Returns
    -------
    T : pandas.DataFrame with columns Name, Description, File — sorted
        alphabetically by Name. Also prints a formatted list to stdout.
    """
    import pandas as pd

    rows = sorted(
        ({"Name": fn.__name__, "Description": _description(fn), "File": fn.__code__.co_filename}
         for fn in _FUNCTIONS),
        key=lambda r: r["Name"],
    )
    T = pd.DataFrame(rows, columns=["Name", "Description", "File"])

    print(f"\n-- G16parser Toolbox -- {len(rows)} function(s) --\n")
    name_w = max(len(r["Name"]) for r in rows)
    for r in rows:
        print(f"  {r['Name']:<{name_w}}  {r['Description']}")
    print()

    return T
