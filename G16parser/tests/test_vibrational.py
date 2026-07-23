import pytest

import G16parser as g16


def test_nmodes_basic(sample_out):
    try:
        nm = g16.g16_nmodes(sample_out)
    except ValueError:
        pytest.skip("sample_out has no vibrational frequency section")
        return

    assert nm.Nmodes > 0
    assert nm.Natoms > 0
    assert nm.disp.shape == (nm.Natoms, 3, nm.Nmodes)
    assert len(nm.freq) == nm.Nmodes


def test_spectra_basic(sample_out):
    try:
        sp = g16.g16_spectra(sample_out)
    except ValueError:
        pytest.skip("sample_out has no vibrational frequency section")
        return

    assert len(sp.freq) > 0
    assert len(sp.x) == len(sp.IR_cont)
    assert (sp.IR_cont >= 0).all()


def test_draw_mode_compatible(sample_out):
    try:
        nm = g16.g16_nmodes(sample_out)
    except ValueError:
        pytest.skip("sample_out has no vibrational frequency section")
        return

    mol = g16.g16_structure(sample_out)
    ax = g16.g16_draw_mode(mol, nm, 1)
    assert ax is not None
