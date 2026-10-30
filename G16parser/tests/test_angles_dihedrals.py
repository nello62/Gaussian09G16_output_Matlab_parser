import re
import warnings

import numpy as np
import pytest

import G16parser as g16
from G16parser._common import Struct

# Ideal geometries with known angles/dihedrals (Angstrom)
SQ = np.array([[0.0, 1.0, 0.0],    # 1
               [0.0, 0.0, 0.0],    # 2
               [1.0, 0.0, 0.0],    # 3
               [1.0, 0.0, 1.0],    # 4 (+90 deg w.r.t. 1)
               [1.0, 0.0, -1.0],   # 5 (-90 deg)
               [1.0, -1.0, 0.0],   # 6 (180 deg, trans)
               [1.0, 1.0, 0.0],    # 7 (0 deg, cis)
               [2.0, 0.0, 0.0]])   # 8 (collinear with 2-3)


def test_angle_known_values():
    _, a = g16.g16_get_angle(SQ, [[1, 2, 3], [2, 3, 8]])
    assert np.allclose(a, [90.0, 180.0])


def test_dihedral_known_values_and_sign():
    _, phi = g16.g16_get_dihedral(SQ, [[1, 2, 3, 4], [1, 2, 3, 5], [1, 2, 3, 6], [1, 2, 3, 7]])
    assert np.allclose(phi, [90.0, -90.0, 180.0, 0.0])
    _, phi_abs = g16.g16_get_dihedral(SQ, [[1, 2, 3, 5]], signed=False)
    assert np.allclose(phi_abs, [90.0])


def test_dihedral_collinear_is_nan():
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        _, phi = g16.g16_get_dihedral(SQ, [2, 3, 8, 4])
    assert np.isnan(phi[0])


def test_table_columns_and_symbols():
    mol = Struct(xyz=SQ, symbols=["H", "C", "C", "H", "H", "H", "H", "C"])
    T, _ = g16.g16_get_angle(mol, [1, 2, 3])
    assert list(T.columns) == ["Atom1", "Sym1", "Atom2", "Sym2", "Atom3", "Sym3", "Angle_deg"]
    assert T.loc[0, "Sym2"] == "C"
    D, _ = g16.g16_get_dihedral(mol, [1, 2, 3, 4])
    assert list(D.columns)[-1] == "Dihedral_deg"


def test_bad_indices_raise():
    with pytest.raises(ValueError):
        g16.g16_get_angle(SQ, [0, 1, 2])
    with pytest.raises(ValueError):
        g16.g16_get_dihedral(SQ, [1, 2, 3, 99])


def test_against_gaussian_optimized_parameters(sample_out):
    """Angles and dihedrals of the final 'Optimized Parameters' table printed
    by Gaussian (if present) are reproduced from the Cartesian geometry."""
    text = open(sample_out, errors="replace").read()
    if "Optimized Parameters" not in text:
        pytest.skip("sample file is not an optimization")
    blk = text[text.rindex("Optimized Parameters"):][:200000]
    A = re.findall(r"! A\d+\s+A\((\d+),(\d+),(\d+)\)\s+(-?\d+\.\d+)", blk)
    D = re.findall(r"! D\d+\s+D\((\d+),(\d+),(\d+),(\d+)\)\s+(-?\d+\.\d+)", blk)
    mol = g16.g16_structure(sample_out)
    if A:
        _, a = g16.g16_get_angle(mol, [[int(x) for x in t[:3]] for t in A])
        assert np.allclose(a, [float(t[3]) for t in A], atol=0.01)
    if D:
        _, phi = g16.g16_get_dihedral(mol, [[int(x) for x in t[:4]] for t in D])
        ref = np.array([float(t[4]) for t in D])
        ok = ~np.isnan(phi)
        diff = (phi[ok] - ref[ok] + 180.0) % 360.0 - 180.0
        assert np.all(np.abs(diff) < 0.01)
