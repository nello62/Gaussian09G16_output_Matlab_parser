import math

import G16parser as g16


def test_energy_basic(sample_out):
    en = g16.g16_energy(sample_out)
    assert en.SCF < 0  # Hartree SCF energies are negative for real molecules
    assert en.filename == sample_out
    if en.has_thermo:
        assert en.E0 > en.SCF  # SCF + (positive) ZPE correction
        assert not math.isnan(en.G_kJ)


def test_dipole_polar_basic(sample_out):
    dp = g16.g16_dipole_polar(sample_out)
    assert dp.mu_tot >= 0
    assert dp.mu_units


def test_orbital_energies_basic(sample_out):
    oe = g16.g16_orbital_energies(sample_out)
    assert oe.HOMO < oe.LUMO
    assert oe.gap > 0
    assert oe.gap_eV > 0
    assert oe.HOMO_eV < oe.LUMO_eV


def test_convergence_runs_without_error(sample_out):
    cv = g16.g16_convergence(sample_out)
    assert cv.Nsteps >= 0
    assert cv.filename == sample_out
    if cv.Nsteps > 0:
        assert len(cv.MaxForce) == cv.Nsteps
