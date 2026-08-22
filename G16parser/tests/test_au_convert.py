import numpy as np
import pytest

import G16parser as g16


@pytest.mark.parametrize("value,quantity,direction,target,expected", [
    (-875.932, "energy", "au2si", "kJ/mol", -2.299759150208e+06),
    (6.3347, "dipole", "si2au", "Debye", 2.492262728538e+00),
    (259.03, "polar", "au2si", "Angstrom3", 3.838428573447e+01),
    (1582.8, "frequency", "si2au", "cm-1", 7.211767433892e-03),
    (0.76, "length", "si2au", "Angstrom", 1.436191855690e+00),
    (1234.5, "energy", "au2si", "cm-1", 2.709414324178e+08),
    (2.5, "beta", "au2si", "esu", 2.159800000000e-32),
    (3.0, "force", "au2si", "nN", 2.471617051152e+02),
    (0.001, "efield", "au2si", "V/Ang", 5.142206751120e-02),
    (298.15, "temp", "si2au", "", 9.414832364531e+07),
    (10.0, "mass", "au2si", "Da", 5.485799090427e-03),
    (1.0, "pressure", "au2si", "GPa", 2.942101575639e+04),
])
def test_au_convert_matches_matlab_reference(value, quantity, direction, target, expected):
    # Reference values from the MATLAB original (AU_convert.m), same
    # CODATA 2022 constants -- agree to ~1e-13 relative difference.
    result = g16.g16_au_convert(value, quantity, direction, target=target, verbose=False)
    assert result == pytest.approx(expected, rel=1e-9)


def test_au_convert_roundtrip():
    au_value = 0.05
    si_value = g16.g16_au_convert(au_value, "energy", "au2si", target="eV", verbose=False)
    back = g16.g16_au_convert(si_value, "energy", "si2au", target="eV", verbose=False)
    assert back == pytest.approx(au_value, rel=1e-12)


def test_au_convert_array_input():
    result = g16.g16_au_convert(np.array([1.0, 2.0, 3.0]), "length", "au2si", target="Angstrom", verbose=False)
    assert result.shape == (3,)
    assert result[1] == pytest.approx(2 * result[0])


def test_au_convert_unknown_quantity_raises():
    with pytest.raises(ValueError):
        g16.g16_au_convert(1.0, "not_a_quantity", "au2si")


def test_au_convert_bad_direction_raises():
    with pytest.raises(ValueError):
        g16.g16_au_convert(1.0, "energy", "sideways")


def test_au_convert_help_returns_none(capsys):
    result = g16.g16_au_convert(None, "help", "")
    assert result is None
    assert "au_convert" in capsys.readouterr().out.lower()


def test_au_table_runs_for_each_section(capsys):
    for section in ("", "energy", "length", "constants", "bogus_section"):
        g16.g16_au_table(section)
    assert "CODATA" in capsys.readouterr().out


@pytest.mark.parametrize("n,unit,expected", [
    (10, "au", 1.0e-3),
    (10, "V/Ang", 5.142206751120e-02),
])
def test_gaussian_field_convert_g2phys(n, unit, expected):
    result = g16.g16_gaussian_field_convert(n, "g2phys", unit=unit, verbose=False)
    assert result == pytest.approx(expected, rel=1e-9)


def test_gaussian_field_convert_phys2g():
    result = g16.g16_gaussian_field_convert(0.0025, "phys2g", unit="au", verbose=False)
    assert result == 25


def test_gaussian_field_convert_roundtrip():
    n = g16.g16_gaussian_field_convert(0.05142206751120, "phys2g", unit="V/Ang", verbose=False)
    assert n == 10


def test_gaussian_field_convert_warns_on_large_rounding():
    with pytest.warns(UserWarning):
        g16.g16_gaussian_field_convert(0.00012345, "phys2g", unit="au", verbose=False)
