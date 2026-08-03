import math
import re
import warnings

import numpy as np

from ._common import Struct, z_to_symbol

_HEADER_RE = re.compile(r'^(.{1,43}?)\s{1,3}(I|R|C)\s+(N=\s*\d+|[-\d.E+]+)')
_SCALAR_TAIL_RE = re.compile(r'(I|R|C)\s+([-\d.E+]+)\s*$')


def _read_lines_raw(filename):
    with open(filename, "r", encoding="utf-8", errors="replace") as f:
        raw = f.read()
    raw = raw.replace("\r\n", "\n").replace("\r", "\n")
    return raw.split("\n")


class _SectionIndex:
    """Indexes every ``.fchk`` section header line, then resolves a
    (possibly fuzzy) keyword to its parsed data -- mirrors the generic
    parser in G09_fchk_read.m/G16_fchk_read.m section for section.
    """

    def __init__(self, lines):
        self.lines = lines
        self.n = len(lines)
        self.names = []
        self.types = []
        self.counts = []
        self.idx = []
        for k, ln in enumerate(lines):
            if len(ln) < 45:
                continue
            m = _HEADER_RE.match(ln)
            if not m:
                continue
            name = m.group(1).strip()
            dtype = m.group(2)
            rest = m.group(3).strip()
            if rest.startswith("N="):
                nvals = int(float(rest[2:].strip()))
            else:
                nvals = 0
            self.names.append(name)
            self.types.append(dtype)
            self.counts.append(nvals)
            self.idx.append(k)

    def read(self, keyword):
        keyword_low = keyword.lower()
        matches = [i for i, name in enumerate(self.names) if name.lower() == keyword_low]
        if not matches:
            matches = [i for i, name in enumerate(self.names) if keyword_low in name.lower()]
        if not matches:
            return None

        mi = matches[-1]
        k0 = self.idx[mi]
        nvals = self.counts[mi]

        if nvals == 0:
            m = _SCALAR_TAIL_RE.search(self.lines[k0])
            if not m:
                return None
            return float(m.group(2))

        vals = []
        k2 = k0 + 1
        while len(vals) < nvals and k2 < self.n:
            ln2 = self.lines[k2].strip()
            if not ln2:
                k2 += 1
                continue
            if len(self.lines[k2]) > 44 and _HEADER_RE.match(self.lines[k2]):
                break
            vals.extend(float(x) for x in ln2.split())
            k2 += 1
        return np.array(vals[:nvals], dtype=float)


def g16_fchk_read(filename, verbose=True):
    """Reads a Gaussian 09/16 formatted checkpoint file (.fchk), generated
    with 'formchk jobname.chk jobname.fchk'. Uses a single generic section
    parser (no hard-coded line counts); despite the g16_ name, this also
    works on G09 .fchk files, since the formatted-checkpoint layout does
    not depend on the Gaussian version.

    Parameters
    ----------
    filename : str
    verbose : bool, default True -- print progress messages while parsing

    Returns
    -------
    data : Struct with (among others) title, method, basis, Nat, charge,
        mult, Nelec, Nalpha, Nbeta, Nbasis, Nbasis_indep, symbols, AN,
        masses, xyz (Angstrom), xyz_bohr, SCF_energy, total_energy,
        virial_ratio, alpha_orb_energies, beta_orb_energies,
        alpha_MO_coeff, mulliken_charges, HOMO_idx, HOMO_eV, LUMO_eV,
        gap_eV, gradient, rms_force, force_const, dipole_au, dipole_D,
        dipole_tot_D, polar_au, polar_iso, polar_aniso, beta_au, beta_vec,
        has_dynamic_beta, dipole_deriv, polar_deriv, filename, and three
        compatibility sub-structs:
            .mol -- compatible with g16_draw_molecule/g16_draw_mode
            .ch  -- compatible with g16_charges_fchk (NOT g16_charges,
                    which only accepts a .log/.out filename)
            .nm  -- compatible with g16_draw_mode(mol, nm, mode_idx)
    """
    with open(filename):
        pass  # raises FileNotFoundError with a clear message if absent
    lines = _read_lines_raw(filename)
    n = len(lines)

    sec = _SectionIndex(lines)

    # -------------------------------------------------------------------
    # Header (lines 1-2)
    # -------------------------------------------------------------------
    title_line = lines[0].strip()
    hdr2 = lines[1].strip()
    tok2 = hdr2.split()
    if len(tok2) >= 3:
        calc_type, method, basis = tok2[0], tok2[1], " ".join(tok2[2:])
    elif len(tok2) == 2:
        method, basis, calc_type = tok2[0], tok2[1], ""
    else:
        method, basis, calc_type = hdr2, "", ""

    # -------------------------------------------------------------------
    # Scalar quantities
    # -------------------------------------------------------------------
    Nat = round(sec.read("Number of atoms"))
    charge = round(sec.read("Charge"))
    mult = round(sec.read("Multiplicity"))
    Nelec = round(sec.read("Number of electrons"))
    Nalpha = round(sec.read("Number of alpha electrons"))
    Nbeta = round(sec.read("Number of beta electrons"))
    Nbasis = round(sec.read("Number of basis functions"))
    Nbasis_indep_raw = sec.read("Number of independent functions")
    Nbasis_indep = round(Nbasis_indep_raw) if Nbasis_indep_raw is not None else None
    SCF_energy = sec.read("SCF Energy")
    total_energy = sec.read("Total Energy")
    virial = sec.read("Virial Ratio")
    rms_force = sec.read("RMS Force")

    if verbose:
        print(f"\n-- g16_fchk_read: {filename} --")
        print(f"  Title  : {title_line}")
        print(f"  Method : {method}  Basis: {basis}")
        print(f"  Nat={Nat}  Charge={charge:+d}  Mult={mult}  Nbasis={Nbasis}")
        print(f"  SCF Energy = {SCF_energy:.10f} Ha")

    # -------------------------------------------------------------------
    # Geometry
    # -------------------------------------------------------------------
    AN_vec = np.round(sec.read("Atomic numbers")).astype(int)
    masses = sec.read("Real atomic weights")
    xyz_bohr_flat = sec.read("Current cartesian coordinates")

    if len(AN_vec) != Nat:
        warnings.warn("g16_fchk_read: Atomic numbers count mismatch.")

    a0_ang = 0.529177210544
    xyz_bohr = xyz_bohr_flat.reshape(Nat, 3)
    xyz_ang = xyz_bohr * a0_ang

    symbols = [z_to_symbol(z) if 1 <= z <= 118 else f"Z{z}" for z in AN_vec]

    if verbose:
        shown = " ".join(symbols[:8])
        more = f" ... ({Nat} total)" if Nat > 8 else ""
        print(f"  Atoms  : {shown}{more}")

    # -------------------------------------------------------------------
    # Electronic structure
    # -------------------------------------------------------------------
    alpha_orb = sec.read("Alpha Orbital Energies")
    beta_orb = sec.read("Beta Orbital Energies")
    alpha_MO = sec.read("Alpha MO coefficients")
    mull_chg = sec.read("Mulliken Charges")

    HOMO_idx = Nalpha
    HOMO_eV = LUMO_eV = gap_eV = math.nan
    if alpha_orb is not None and len(alpha_orb) >= HOMO_idx + 1:
        Eh2eV = 27.211386
        HOMO_eV = alpha_orb[HOMO_idx - 1] * Eh2eV
        LUMO_eV = alpha_orb[HOMO_idx] * Eh2eV
        gap_eV = LUMO_eV - HOMO_eV
        if verbose:
            print(f"  HOMO   = {HOMO_eV:.4f} eV   LUMO = {LUMO_eV:.4f} eV   Gap = {gap_eV:.4f} eV")

    # -------------------------------------------------------------------
    # Forces
    # -------------------------------------------------------------------
    N3 = 3 * Nat
    grad_raw = sec.read("Cartesian Gradient")
    gradient = np.full((Nat, 3), np.nan)
    if grad_raw is not None and len(grad_raw) == N3:
        gradient = grad_raw.reshape(Nat, 3)

    fc_raw = sec.read("Cartesian Force Constants")
    force_const = np.full((N3, N3), np.nan)
    if fc_raw is not None and len(fc_raw) == N3 * (N3 + 1) // 2:
        k_fc = 0
        for row in range(N3):
            for col in range(row + 1):
                force_const[row, col] = fc_raw[k_fc]
                force_const[col, row] = fc_raw[k_fc]
                k_fc += 1

    # -------------------------------------------------------------------
    # Dipole moment
    # -------------------------------------------------------------------
    dip_raw = sec.read("Dipole Moment")
    dip_au = np.array([math.nan] * 3)
    dip_D = np.array([math.nan] * 3)
    dip_tot_D = math.nan
    if dip_raw is not None and len(dip_raw) == 3:
        au2D = 2.541747
        dip_au = dip_raw
        dip_D = dip_au * au2D
        dip_tot_D = float(np.linalg.norm(dip_D))
        if verbose:
            print(f"  Dipole  = ({dip_D[0]:.4f}, {dip_D[1]:.4f}, {dip_D[2]:.4f}) D   "
                  f"|mu| = {dip_tot_D:.4f} D")

    # -------------------------------------------------------------------
    # Polarisability
    # -------------------------------------------------------------------
    pol_raw = sec.read("Polarizability")
    pol_au = np.full((3, 3), np.nan)
    pol_iso = pol_aniso = math.nan
    if pol_raw is not None and len(pol_raw) == 6:
        v = pol_raw
        pol_au = np.array([[v[0], v[1], v[3]],
                            [v[1], v[2], v[4]],
                            [v[3], v[4], v[5]]])
        pol_iso = float(np.trace(pol_au) / 3)
        pol_aniso = float(math.sqrt(0.5 * (
            (pol_au[0, 0] - pol_au[1, 1]) ** 2 +
            (pol_au[1, 1] - pol_au[2, 2]) ** 2 +
            (pol_au[2, 2] - pol_au[0, 0]) ** 2 +
            6 * (pol_au[0, 1] ** 2 + pol_au[0, 2] ** 2 + pol_au[1, 2] ** 2))))
        if verbose:
            print(f"  alpha_iso   = {pol_iso:.4f} au   alpha_aniso = {pol_aniso:.4f} au")

    # -------------------------------------------------------------------
    # First hyperpolarisability
    # -------------------------------------------------------------------
    hyper_raw = sec.read("HyperPolarizability")
    beta_fields = ["xxx", "xxy", "xyy", "yyy", "xxz", "xyz", "yyz", "xzz", "yzz", "zzz"]
    beta = {f: math.nan for f in beta_fields}
    beta_vec = math.nan
    if hyper_raw is not None and len(hyper_raw) >= 10:
        for f, val in zip(beta_fields, hyper_raw[:10]):
            beta[f] = val
        bx = beta["xxx"] + beta["xyy"] + beta["xzz"]
        by = beta["xxy"] + beta["yyy"] + beta["yzz"]
        bz = beta["xxz"] + beta["yyz"] + beta["zzz"]
        beta_vec = math.sqrt(bx ** 2 + by ** 2 + bz ** 2)
        if verbose:
            print(f"  |beta_vec| = {beta_vec:.2f} au")

    has_dynamic_beta = any(re.match(r"^Beta\(.*\)$", name) for name in sec.names)
    if has_dynamic_beta and math.isnan(beta_vec):
        warnings.warn(
            f"{filename} contains a frequency-dependent hyperpolarisability section "
            "(dynamic Beta, e.g. from CPHF=RdFreq) that this parser does not yet "
            "decode -- beta_au/beta_vec are NaN even though hyperpolarisability "
            "data is present. See data.has_dynamic_beta."
        )

    # -------------------------------------------------------------------
    # Dipole / polarisability derivatives
    # -------------------------------------------------------------------
    dip_deriv_raw = sec.read("Dipole Derivatives")
    dip_deriv = np.full((3, N3), np.nan)
    if dip_deriv_raw is not None and len(dip_deriv_raw) == 3 * N3:
        dip_deriv = dip_deriv_raw.reshape(N3, 3).T

    pol_deriv_raw = sec.read("Polarizability Derivatives")
    pol_deriv = np.full((6, N3), np.nan)
    if pol_deriv_raw is not None and len(pol_deriv_raw) == 6 * N3:
        pol_deriv = pol_deriv_raw.reshape(N3, 6).T

    # -------------------------------------------------------------------
    # Assemble output struct
    # -------------------------------------------------------------------
    data = Struct(
        title=title_line, method=method, basis=basis, calc_type=calc_type,
        Nat=Nat, charge=charge, mult=mult, Nelec=Nelec, Nalpha=Nalpha, Nbeta=Nbeta,
        Nbasis=Nbasis, Nbasis_indep=Nbasis_indep,
        symbols=symbols, AN=AN_vec, masses=masses, xyz=xyz_ang, xyz_bohr=xyz_bohr,
        SCF_energy=SCF_energy, total_energy=total_energy, virial_ratio=virial,
        rms_force=rms_force,
        alpha_orb_energies=alpha_orb, beta_orb_energies=beta_orb,
        alpha_MO_coeff=alpha_MO, mulliken_charges=mull_chg,
        HOMO_idx=HOMO_idx, HOMO_eV=HOMO_eV, LUMO_eV=LUMO_eV, gap_eV=gap_eV,
        gradient=gradient, force_const=force_const,
        dipole_au=dip_au, dipole_D=dip_D, dipole_tot_D=dip_tot_D,
        polar_au=pol_au, polar_iso=pol_iso, polar_aniso=pol_aniso,
        beta_au=Struct(**beta), beta_vec=beta_vec, has_dynamic_beta=has_dynamic_beta,
        dipole_deriv=dip_deriv, polar_deriv=pol_deriv,
        filename=filename,
    )

    # -------------------------------------------------------------------
    # .mol -- compatible with g16_draw_molecule / g16_draw_mode
    # -------------------------------------------------------------------
    data.mol = Struct(
        symbols=symbols, xyz=xyz_ang, Z=AN_vec, Natoms=Nat,
        step=1, n_steps=1, orientation="fchk (Input orientation)",
        filename=filename,
    )

    # -------------------------------------------------------------------
    # .ch -- compatible with g16_charges_fchk (NOT g16_charges, which only
    #        accepts a .log/.out filename, not a pre-parsed struct)
    # -------------------------------------------------------------------
    data.ch = Struct(
        symbols=symbols, charges=mull_chg, charges_H=None,
        sum_q=float(np.sum(mull_chg)) if mull_chg is not None else math.nan,
        type="Mulliken", label="Mulliken Charges (from .fchk)",
        Natoms=Nat, filename=filename,
    )

    # -------------------------------------------------------------------
    # .nm -- compatible with g16_draw_mode(mol, nm, mode_idx)
    # -------------------------------------------------------------------
    nm = Struct(
        Nmodes=0, Natoms=Nat, has_Raman=False,
        freq=np.array([]), IR=np.array([]), Raman=np.array([]),
        disp=np.zeros((Nat, 3, 0)), symmetry=[], redmass=None, frcconst=None,
        filename=filename,
    )

    if not np.any(np.isnan(force_const)) and masses is not None and len(masses) == Nat:
        m_vec = np.repeat(masses, 3)   # [3Nat] amu

        M_tot = np.sum(masses)
        com = (masses @ xyz_bohr) / M_tot
        xyz_c = xyz_bohr - com

        D_tr = np.zeros((N3, 6))
        for ii in range(Nat):
            sqm = math.sqrt(masses[ii])
            xi, yi, zi = xyz_c[ii]
            D_tr[3 * ii, 0] = sqm
            D_tr[3 * ii + 1, 1] = sqm
            D_tr[3 * ii + 2, 2] = sqm
            D_tr[3 * ii, 3] = 0
            D_tr[3 * ii + 1, 3] = -zi * sqm
            D_tr[3 * ii + 2, 3] = yi * sqm
            D_tr[3 * ii, 4] = zi * sqm
            D_tr[3 * ii + 1, 4] = 0
            D_tr[3 * ii + 2, 4] = -xi * sqm
            D_tr[3 * ii, 5] = -yi * sqm
            D_tr[3 * ii + 1, 5] = xi * sqm
            D_tr[3 * ii + 2, 5] = 0

        Q_tr, R_tr = np.linalg.qr(D_tr)
        tr_keep = np.abs(np.diag(R_tr)) > 1e-6
        D_orth = Q_tr[:, tr_keep]
        n_tr = D_orth.shape[1]

        P_vib = np.eye(N3) - D_orth @ D_orth.T

        mw = 1.0 / np.sqrt(m_vec)
        Fmw = force_const * np.outer(mw, mw)
        Fmw = (Fmw + Fmw.T) / 2

        Fmw_proj = P_vib @ Fmw @ P_vib
        Fmw_proj = (Fmw_proj + Fmw_proj.T) / 2

        lam, V = np.linalg.eigh(Fmw_proj)

        Eh_J = 4.3597447222060e-18
        a0_m = 5.29177210544e-11
        amu_kg = 1.66053906892e-27
        c_cms = 2.99792458e10
        hess2cm1 = math.sqrt(Eh_J / (a0_m ** 2 * amu_kg)) / (2 * math.pi * c_cms)

        freq_cm1 = np.sign(lam) * np.sqrt(np.abs(lam)) * hess2cm1

        sort_idx = np.argsort(freq_cm1)
        freq_sorted = freq_cm1[sort_idx]

        L_cart = mw[:, None] * V
        col_norms = np.sqrt(np.sum(L_cart ** 2, axis=0))
        col_norms[col_norms == 0] = 1
        L_cart = L_cart / col_norms
        L_sorted = L_cart[:, sort_idx]

        freq_vib = freq_sorted[n_tr:]
        L_vib_cn = L_sorted[:, n_tr:]
        Nmodes_vib = N3 - n_tr

        L_mn_sorted = mw[:, None] * V[:, sort_idx]
        L_vib_mn = L_mn_sorted[:, n_tr:]

        disp_all = L_vib_cn.reshape(3, Nat, Nmodes_vib, order="F")
        disp_all = np.transpose(disp_all, (1, 0, 2))

        # dip_deriv/pol_deriv are already (3,N3)/(6,N3): row = Cartesian
        # dipole (x/y/z) or alpha (xx,xy,yy,xz,yz,zz) component, column =
        # displacement direction -- used directly, with no further
        # reshape (the equivalent MATLAB code once re-reshaped these
        # already-correctly-shaped arrays to (N3,3)/(N3,6) and back,
        # which does NOT undo itself -- it reinterprets the same
        # column-major buffer under a different shape, scrambling
        # components across displacements. Caught by comparing against
        # the same molecule's IR/Raman as printed directly by Gaussian
        # in the corresponding .out file; fixed in G09/G16_fchk_read.m).
        C_IR = 42.2561 * (2.541747 / 0.529177) ** 2
        IR_int = np.zeros(Nmodes_vib)
        if not np.any(np.isnan(dip_deriv)):
            dmu_dQ = dip_deriv @ L_vib_mn
            IR_int = np.sum(dmu_dQ ** 2, axis=0) * C_IR

        C_Ra = (0.148185 / 0.529177) ** 2
        Raman_int = np.zeros(Nmodes_vib)
        has_raman = False
        if not np.any(np.isnan(pol_deriv)):
            da_dQ = pol_deriv @ L_vib_mn
            d_iso = (da_dQ[0] + da_dQ[2] + da_dQ[5]) / 3
            d_an2 = ((da_dQ[0] - da_dQ[2]) ** 2 +
                     (da_dQ[2] - da_dQ[5]) ** 2 +
                     (da_dQ[5] - da_dQ[0]) ** 2) / 2 + \
                    3 * (da_dQ[1] ** 2 + da_dQ[3] ** 2 + da_dQ[4] ** 2)
            Raman_int = (45 * d_iso ** 2 + 7 * d_an2) * C_Ra
            has_raman = True

        nm = Struct(
            Nmodes=Nmodes_vib, Natoms=Nat, has_Raman=has_raman,
            freq=freq_vib, IR=IR_int, Raman=Raman_int, disp=disp_all,
            symmetry=["?"] * Nmodes_vib, redmass=None, frcconst=None,
            filename=filename,
        )

        if verbose:
            print(f"  Normal modes: {Nmodes_vib} vibrational  ({n_tr} TR projected out)")
            real_mask = freq_vib > 10
            if np.any(real_mask):
                first_real = np.argmax(real_mask)
                print(f"  Lowest freq  : {freq_vib[first_real]:.1f} cm-1 (mode {first_real + 1})")
                print(f"  Highest freq : {freq_vib[-1]:.1f} cm-1")
            n_imag = int(np.sum(freq_vib < -5))
            if n_imag > 0:
                print(f"  Imaginary (<-5 cm-1): {n_imag}  (saddle point / TS geometry)")

    data.nm = nm

    if verbose:
        print("\n  Compatibility sub-structs ready:")
        print("    data.mol  -> g16_draw_molecule(data.mol)")
        print("    data.ch   -> g16_charges_fchk(data.mol, data.ch)")
        print("    data.nm   -> g16_draw_mode(data.mol, data.nm, mode_idx)")
        print("  Done.\n")

    return data
