import numpy as np

from ._common import Struct, read_lines, extract_floats, HARTREE_TO_EV


def _push_block(blocks, a_occ, a_virt, b_occ, b_virt):
    if a_occ or a_virt:
        blocks.append({
            "alpha_occ": np.array(a_occ), "alpha_virt": np.array(a_virt),
            "beta_occ": np.array(b_occ), "beta_virt": np.array(b_virt),
        })


def g16_orbital_energies(filename, step="last", lines=None):
    """Extracts molecular orbital energies (HOMO/LUMO and the full
    occupied/virtual spectrum) from a Gaussian 16 output file.

    Reads the "Alpha  occ. eigenvalues --" / "Alpha virt. eigenvalues --"
    (and, for open-shell calculations, the matching "Beta" lines) printed
    by Gaussian's population analysis. These blocks repeat once per SCF
    calculation in the file (e.g. once per optimisation step); `step`
    selects which block to report on, exactly like g16_energy/g16_structure.

    Returns
    -------
    oe : Struct — alpha_occ, alpha_virt, beta_occ, beta_virt (np.ndarray),
        has_beta, HOMO, LUMO, gap, HOMO_alpha, LUMO_alpha, HOMO_beta,
        LUMO_beta, HOMO_eV, LUMO_eV, gap_eV, step, Nsteps, filename.
    """
    if lines is None:
        lines = read_lines(filename)

    blocks = []
    cur = {"alpha_occ": [], "alpha_virt": [], "beta_occ": [], "beta_virt": []}
    last_kind = ""

    for ln in lines:
        if "Alpha  occ. eigenvalues" in ln or "Alpha occ. eigenvalues" in ln:
            if last_kind in ("alpha_virt", "beta_virt"):
                _push_block(blocks, cur["alpha_occ"], cur["alpha_virt"],
                            cur["beta_occ"], cur["beta_virt"])
                cur = {"alpha_occ": [], "alpha_virt": [], "beta_occ": [], "beta_virt": []}
            cur["alpha_occ"].extend(extract_floats(ln))
            last_kind = "alpha_occ"
        elif "Alpha virt. eigenvalues" in ln:
            cur["alpha_virt"].extend(extract_floats(ln))
            last_kind = "alpha_virt"
        elif "Beta  occ. eigenvalues" in ln or "Beta occ. eigenvalues" in ln:
            cur["beta_occ"].extend(extract_floats(ln))
            last_kind = "beta_occ"
        elif "Beta virt. eigenvalues" in ln:
            cur["beta_virt"].extend(extract_floats(ln))
            last_kind = "beta_virt"

    _push_block(blocks, cur["alpha_occ"], cur["alpha_virt"], cur["beta_occ"], cur["beta_virt"])

    if not blocks:
        raise ValueError(f"g16_orbital_energies: no orbital eigenvalue block found in {filename}")

    nsteps = len(blocks)
    if isinstance(step, str):
        if step.lower() == "last":
            si = nsteps
        elif step.lower() == "first":
            si = 1
        else:
            raise ValueError("g16_orbital_energies: step must be 'first', 'last', or an integer.")
    else:
        si = round(step)
        if si < 1 or si > nsteps:
            raise ValueError(f"g16_orbital_energies: step {si} out of range [1, {nsteps}].")

    blk = blocks[si - 1]
    if blk["alpha_occ"].size == 0 or blk["alpha_virt"].size == 0:
        raise ValueError(f"g16_orbital_energies: incomplete eigenvalue block (step {si}) in {filename}")

    HOMO_a = blk["alpha_occ"][-1]
    LUMO_a = blk["alpha_virt"][0]
    has_beta = blk["beta_occ"].size > 0 and blk["beta_virt"].size > 0

    if has_beta:
        HOMO_b = blk["beta_occ"][-1]
        LUMO_b = blk["beta_virt"][0]
        HOMO = max(HOMO_a, HOMO_b)
        LUMO = min(LUMO_a, LUMO_b)
    else:
        HOMO_b = LUMO_b = None
        HOMO, LUMO = HOMO_a, LUMO_a

    gap = LUMO - HOMO

    oe = Struct(
        alpha_occ=blk["alpha_occ"], alpha_virt=blk["alpha_virt"],
        beta_occ=blk["beta_occ"], beta_virt=blk["beta_virt"], has_beta=has_beta,
        HOMO=HOMO, LUMO=LUMO, gap=gap,
        HOMO_alpha=HOMO_a, LUMO_alpha=LUMO_a, HOMO_beta=HOMO_b, LUMO_beta=LUMO_b,
        HOMO_eV=HOMO * HARTREE_TO_EV, LUMO_eV=LUMO * HARTREE_TO_EV, gap_eV=gap * HARTREE_TO_EV,
        step=si, Nsteps=nsteps, filename=filename,
    )

    print(f"\n── g16_orbital_energies (step {si}/{nsteps}): {filename} ──")
    print(f"  HOMO = {HOMO:+.6f} Ha  ({oe.HOMO_eV:+.4f} eV)")
    print(f"  LUMO = {LUMO:+.6f} Ha  ({oe.LUMO_eV:+.4f} eV)")
    print(f"  Gap  =  {gap:.6f} Ha  ( {oe.gap_eV:.4f} eV)")
    if has_beta:
        print(f"  (open-shell: alpha HOMO/LUMO = {HOMO_a:+.6f} / {LUMO_a:+.6f} Ha, "
              f"beta HOMO/LUMO = {HOMO_b:+.6f} / {LUMO_b:+.6f} Ha)")
    print()

    return oe
