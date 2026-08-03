import numpy as np

import G16parser as g16


def test_fchk_read_basic(sample_fchk):
    data = g16.g16_fchk_read(sample_fchk, verbose=False)
    assert data.Nat > 0
    assert len(data.symbols) == data.Nat
    assert data.xyz.shape == (data.Nat, 3)
    assert data.filename == sample_fchk


def test_fchk_read_compat_structs(sample_fchk):
    data = g16.g16_fchk_read(sample_fchk, verbose=False)

    assert data.mol.Natoms == data.Nat
    assert data.mol.xyz.shape == (data.Nat, 3)

    assert data.ch.Natoms == data.Nat
    assert data.ch.type == "Mulliken"

    assert data.nm.Natoms == data.Nat
    if data.nm.Nmodes > 0:
        assert data.nm.disp.shape == (data.Nat, 3, data.nm.Nmodes)
        assert len(data.nm.freq) == data.nm.Nmodes
        assert len(data.nm.IR) == data.nm.Nmodes


def test_fchk_read_draw_molecule_compatible(sample_fchk):
    data = g16.g16_fchk_read(sample_fchk, verbose=False)
    ax = g16.g16_draw_molecule(data.mol)
    assert ax is not None


def test_charges_fchk_basic(sample_fchk):
    data = g16.g16_fchk_read(sample_fchk, verbose=False)
    ch = g16.g16_charges_fchk(data.mol, data.ch, plot=False)
    assert ch.Natoms == data.Nat
    assert len(ch.charges) == data.Nat


def test_charges_fchk_missing_field_raises():
    from G16parser._common import Struct
    bad_mol = Struct(symbols=["H"], xyz=np.zeros((1, 3)))  # no Natoms
    ch = Struct(charges=[0.0], symbols=["H"], type="Mulliken", Natoms=1)
    try:
        g16.g16_charges_fchk(bad_mol, ch, plot=False)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_fchk_ir_raman_matches_out_ground_truth(sample_fchk_pair):
    # Regression test for a real bug found while porting: G09_fchk_read.m/
    # G16_fchk_read.m used to reshape the already-correctly-shaped
    # dip_deriv/pol_deriv arrays a second time before the IR/Raman
    # intensity calculation, which silently scrambled the data (still
    # produced plausible-looking but wrong numbers). Caught by comparing
    # against the IR/Raman Gaussian itself prints in the .out file for
    # the same calculation. This Python port never had the bug (fixed
    # directly), but this test guards against it ever being (re)introduced.
    fchk_path, out_path = sample_fchk_pair

    data = g16.g16_fchk_read(fchk_path, verbose=False)
    nm_out = g16.g16_nmodes(out_path)

    assert data.nm.Nmodes == nm_out.Nmodes
    np.testing.assert_allclose(data.nm.freq, nm_out.freq, atol=0.1)
    np.testing.assert_allclose(data.nm.IR, nm_out.IR, rtol=0.01, atol=0.01)
    if nm_out.has_Raman and data.nm.has_Raman:
        np.testing.assert_allclose(data.nm.Raman, nm_out.Raman, rtol=0.01, atol=0.01)


def test_spectra_nm_matches_spectra_ground_truth(sample_fchk_pair):
    # g16_spectra_nm builds the same Lorentzian-broadened continuum as
    # g16_spectra, but from an already-parsed nm struct (e.g. data.nm from
    # g16_fchk_read) instead of re-reading a .log/.out file. Cross-check
    # against g16_spectra called directly on the matching .out file.
    fchk_path, out_path = sample_fchk_pair

    sp_ref = g16.g16_spectra(out_path, FWHM=12, xmin=0, xmax=4000, dx=2)

    data = g16.g16_fchk_read(fchk_path, verbose=False)
    sp_nm = g16.g16_spectra_nm(data.nm, FWHM=12, xmin=0, xmax=4000, dx=2)

    assert sp_ref.Nmodes == sp_nm.Nmodes
    assert sp_nm.has_Raman == sp_ref.has_Raman
    np.testing.assert_array_equal(sp_ref.x, sp_nm.x)
    np.testing.assert_allclose(sp_ref.IR_cont, sp_nm.IR_cont, rtol=0.01, atol=0.01)
    if sp_ref.has_Raman:
        np.testing.assert_allclose(sp_ref.Raman_cont, sp_nm.Raman_cont, rtol=0.01, atol=0.01)

    # Also works directly with g16_nmodes's output (not just fchk_read)
    nm2 = g16.g16_nmodes(out_path)
    sp_nm2 = g16.g16_spectra_nm(nm2, FWHM=12, xmin=0, xmax=4000, dx=2)
    np.testing.assert_allclose(sp_ref.IR_cont, sp_nm2.IR_cont, atol=1e-9)


def test_spectra_nm_normalize(sample_fchk):
    data = g16.g16_fchk_read(sample_fchk, verbose=False)
    sp = g16.g16_spectra_nm(data.nm, normalize=True)
    assert abs(sp.IR_cont.max() - 1) < 1e-9
    if sp.has_Raman:
        assert abs(sp.Raman_cont.max() - 1) < 1e-9


def test_spectra_nm_missing_field_raises():
    from G16parser._common import Struct
    try:
        g16.g16_spectra_nm(Struct(foo=1))
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_spectra_nm_zero_modes_raises():
    from G16parser._common import Struct
    try:
        g16.g16_spectra_nm(Struct(freq=np.array([]), IR=np.array([]), Nmodes=0))
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_spectra_nm_no_raman():
    from G16parser._common import Struct
    nm = Struct(freq=np.array([100.0, 200.0]), IR=np.array([1.0, 2.0]),
                Nmodes=2, filename="test.out")
    sp = g16.g16_spectra_nm(nm)
    assert sp.has_Raman is False
    assert len(sp.Raman) == 0
    assert len(sp.Raman_cont) == 0
