import os

import numpy as np
import pandas as pd

from .read_all import g16_read_all
from .orbital_energies import g16_orbital_energies
from .write_report import g16_write_report

_EXTENSIONS = (".log", ".out")


def g16_batch_read_all(folder, recursive=False, write_reports=False, save_as=""):
    """Runs g16_read_all over every Gaussian 16 output file in a folder and
    aggregates the key results into one summary DataFrame.

    Port of G16_batch_read_all.m.

    Scans `folder` for .log/.out files (case-insensitive), runs
    g16_read_all (plus g16_orbital_energies, for the HOMO/LUMO gap) on
    each, and returns one row per file. A file that fails to parse (e.g.
    an incomplete/crashed job, or a Gaussian 09 file run through this
    Gaussian-16-only package by mistake) does not stop the batch: its row
    is filled with NaN and the error message is recorded in the Status
    column instead.

    Parameters
    ----------
    folder : str
    recursive : bool, optional — also scan subfolders (default False)
    write_reports : bool, optional — write a g16_write_report .txt next to
        each source file (default False)
    save_as : str, optional — path to save the summary table (.csv or
        .xlsx, inferred from the extension)

    Returns
    -------
    T : pandas.DataFrame with one row per file and columns: File, Natoms,
        SCF_Hartree, E0_Hartree, G_kJmol, mu_tot, mu_units, HOMO_eV,
        LUMO_eV, Gap_eV, Status
    """
    if not os.path.isdir(folder):
        raise ValueError(f"g16_batch_read_all: folder not found: {folder}")

    files = []
    if recursive:
        for root, _dirs, names in os.walk(folder):
            for name in names:
                if name.lower().endswith(_EXTENSIONS):
                    files.append(os.path.join(root, name))
    else:
        for name in sorted(os.listdir(folder)):
            full = os.path.join(folder, name)
            if os.path.isfile(full) and name.lower().endswith(_EXTENSIONS):
                files.append(full)

    if not files:
        raise ValueError(f"g16_batch_read_all: no .log/.out files found in {folder}")

    rows = []
    n_ok = 0
    for fullpath in files:
        row = {
            "File": fullpath, "Natoms": np.nan, "SCF_Hartree": np.nan,
            "E0_Hartree": np.nan, "G_kJmol": np.nan, "mu_tot": np.nan,
            "mu_units": "", "HOMO_eV": np.nan, "LUMO_eV": np.nan,
            "Gap_eV": np.nan, "Status": "",
        }
        try:
            Tk = g16_read_all(fullpath)
            oe = g16_orbital_energies(fullpath)

            row["Natoms"] = Tk.structure.Natoms
            row["SCF_Hartree"] = Tk.energy.SCF
            row["E0_Hartree"] = Tk.energy.E0
            row["G_kJmol"] = Tk.energy.G_kJ
            row["mu_tot"] = Tk.dipolar.mu_tot
            row["mu_units"] = Tk.dipolar.mu_units
            row["HOMO_eV"] = oe.HOMO_eV
            row["LUMO_eV"] = oe.LUMO_eV
            row["Gap_eV"] = oe.gap_eV
            row["Status"] = "ok"
            n_ok += 1

            if write_reports:
                fdir, fname = os.path.split(fullpath)
                base, _ = os.path.splitext(fname)
                g16_write_report(Tk, os.path.join(fdir, f"{base}_report.txt"))
        except Exception as err:
            row["Status"] = str(err)

        rows.append(row)

    T = pd.DataFrame(rows, columns=[
        "File", "Natoms", "SCF_Hartree", "E0_Hartree", "G_kJmol", "mu_tot",
        "mu_units", "HOMO_eV", "LUMO_eV", "Gap_eV", "Status",
    ])

    if save_as:
        if save_as.lower().endswith(".xlsx"):
            T.to_excel(save_as, index=False)
        else:
            T.to_csv(save_as, index=False)

    print("\n── g16_batch_read_all ──")
    print(f"  Folder     : {folder}")
    print(f"  Files found: {len(files)}")
    print(f"  Succeeded  : {n_ok}")
    print(f"  Failed     : {len(files) - n_ok}")
    if save_as:
        print(f"  Saved to   : {save_as}")
    print()

    return T
