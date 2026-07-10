import datetime
import os


def g16_write_report(T, outfile=None):
    """Writes a human-readable text report from a g16_read_all Struct.

    Python port of G16_write_report.m: writes a formatted summary of every
    field in T (the Struct returned by g16_read_all) to a .txt file, named
    after the source Gaussian file by default (e.g. 'zeatin.out' ->
    'zeatin_report.txt' in the current folder).

    Sections included (only if the corresponding attribute is present on
    T): route, charge/multiplicity, molecular structure, energetics,
    dipole moment and polarisability, atomic charges, vibrational modes,
    and a summary of the simulated IR/Raman spectra (the full continuum
    arrays are not dumped — use T.spectra.x / .IR_cont / .Raman_cont
    directly for that).

    Parameters
    ----------
    T : Struct — from g16_read_all
    outfile : str, optional — output path (default: "<source>_report.txt")

    Returns
    -------
    outfile : str — the path written to.
    """
    structure = getattr(T, "structure", None)

    if not outfile:
        src = getattr(structure, "filename", None) if structure is not None else None
        if not src:
            outfile = "G16_report.txt"
        else:
            name = os.path.splitext(os.path.basename(src))[0]
            outfile = f"{name}_report.txt"

    out = []

    def w(line=""):
        out.append(line)

    w("=" * 64)
    w("  Gaussian 16 Calculation Report")
    w("=" * 64)
    if structure is not None and getattr(structure, "filename", None):
        w(f"Source file : {structure.filename}")
    w(f"Generated   : {datetime.datetime.now():%d-%b-%Y %H:%M:%S}")
    w()

    route = getattr(T, "route", None)
    if route is not None:
        w("--- Route section " + "-" * 46)
        w(route)
        w()

    chargemol = getattr(T, "chargemol", None)
    if chargemol is not None:
        w("--- Charge / multiplicity " + "-" * 38)
        w(f"Total charge       : {chargemol.charge}")
        w(f"Spin multiplicity  : {chargemol.mol}")
        w()

    if structure is not None:
        s = structure
        w("--- Molecular structure " + "-" * 40)
        w(f"Atoms       : {s.Natoms}")
        if getattr(s, "orientation", None):
            w(f"Orientation : {s.orientation}")
        if getattr(s, "step", None) is not None:
            w(f"Step        : {s.step}")
        w()
        w(f"{'Idx':<6} {'Sym':<4} {'X (A)':>12} {'Y (A)':>12} {'Z (A)':>12}")
        for i in range(s.Natoms):
            x, y, z = s.xyz[i]
            w(f"{i + 1:<6} {s.symbols[i]:<4} {x:12.6f} {y:12.6f} {z:12.6f}")
        w()

    energy = getattr(T, "energy", None)
    if energy is not None:
        e = energy
        w("--- Energetics " + "-" * 49)
        w(f"Method          : {e.method}")
        w(f"SCF energy      : {e.SCF:.8f} Hartree")
        if e.has_thermo:
            w(f"ZPE correction  : {e.ZPE_corr:.8f} Hartree  ({e.ZPE_kJ:.3f} kJ/mol)")
            w(f"Thermal U corr. : {e.U_corr:.8f} Hartree")
            w(f"Thermal H corr. : {e.H_corr:.8f} Hartree")
            w(f"Thermal G corr. : {e.G_corr:.8f} Hartree")
            w(f"E0 (SCF+ZPE)    : {e.E0:.8f} Hartree")
            w(f"U               : {e.U:.8f} Hartree")
            w(f"H               : {e.H:.8f} Hartree")
            w(f"G               : {e.G:.8f} Hartree")
            w(f"T, P            : {e.T:.2f} K, {e.P:.4f} atm")
        else:
            w("(no thermochemistry data - opt/single-point job without freq)")
        w()

    dipolar = getattr(T, "dipolar", None)
    if dipolar is not None:
        d = dipolar
        w("--- Dipole moment and polarisability " + "-" * 27)
        w(f"Dipole (mu_x, mu_y, mu_z) : {d.mu_x:.6f}  {d.mu_y:.6f}  {d.mu_z:.6f}  [{d.mu_units}]")
        w(f"Dipole magnitude          : {d.mu_tot:.6f} {d.mu_units}")
        w(f"Alpha isotropic           : {d.alpha_iso:.6f} {d.alpha_units}")
        w(f"Alpha anisotropy          : {d.alpha_aniso:.6f} {d.alpha_units}")
        if getattr(d, "N_dyn", 0) > 0:
            w()
            w("Dynamic polarisability Alpha(-w;w):")
            w(f"{'Lambda (nm)':<14} {'Freq (au)':<14} {'Iso':>14} {'Aniso':>14}")
            for ad in d.alpha_dyn:
                w(f"{ad.lambda_nm:<14.2f} {ad.freq_au:<14.6f} {ad.iso:14.6f} {ad.aniso:14.6f}")
        w()

    charge = getattr(T, "charge", None)
    if charge is not None:
        c = charge
        w(f"--- Atomic charges ({c.type}) " + "-" * 30)
        w(f"{'Idx':<6} {'Sym':<4} {'Charge (e)':>14}")
        for i in range(c.Natoms):
            w(f"{i + 1:<6} {c.symbols[i]:<4} {c.charges[i]:14.6f}")
        w(f"Sum of charges : {c.sum_q:.6f}")
        if getattr(c, "dipole", None) is not None:
            mu = c.dipole
            w(f"Dipole (from charges overlay) : {mu[0]:.6f}  {mu[1]:.6f}  {mu[2]:.6f}  Debye")
        w()

    nmodes = getattr(T, "nmodes", None)
    if nmodes is not None:
        nm = nmodes
        w("--- Vibrational normal modes " + "-" * 35)
        w(f"Number of modes : {nm.Nmodes}")
        w()
        if nm.has_Raman:
            w(f"{'Mode':<6} {'Freq(cm-1)':>12} {'IR':>10} {'Raman':>12} {'RedMass':>10} {'Sym':>8}")
        else:
            w(f"{'Mode':<6} {'Freq(cm-1)':>12} {'IR':>10} {'RedMass':>10} {'Sym':>8}")
        for k in range(nm.Nmodes):
            sym = nm.symmetry[k] if k < len(nm.symmetry) else ""
            if nm.has_Raman:
                w(f"{k + 1:<6} {nm.freq[k]:12.2f} {nm.IR[k]:10.2f} {nm.Raman[k]:12.2f} "
                  f"{nm.redmass[k]:10.4f} {sym:>8}")
            else:
                w(f"{k + 1:<6} {nm.freq[k]:12.2f} {nm.IR[k]:10.2f} {nm.redmass[k]:10.4f} {sym:>8}")
        w()

    spectra = getattr(T, "spectra", None)
    if spectra is not None:
        sp = spectra
        w("--- Simulated IR/Raman spectra " + "-" * 33)
        w(f"FWHM used     : {sp.FWHM:.2f} cm^-1")
        w(f"Grid range    : {sp.x[0]:.1f} - {sp.x[-1]:.1f} cm^-1 ({len(sp.x)} points)")
        w(f"Raman present : {sp.has_Raman}")
        w("(full continuum arrays not dumped here - see T.spectra.x / .IR_cont / .Raman_cont)")
        w()

    with open(outfile, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")

    print(f"g16_write_report: report written to {outfile}")
    return outfile
