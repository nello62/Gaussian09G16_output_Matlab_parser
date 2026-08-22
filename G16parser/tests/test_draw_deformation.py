import copy
import warnings

import numpy as np
import pytest

import G16parser as g16


def test_draw_deformation_basic(sample_out):
    mol1 = g16.g16_structure(sample_out)
    mol2 = copy.deepcopy(mol1)
    rng = np.random.default_rng(0)
    mol2.xyz = mol1.xyz + rng.normal(scale=0.02, size=mol1.xyz.shape)

    ax = g16.g16_draw_deformation(mol1, mol2, scale=20)
    assert ax is not None


def test_draw_deformation_overlay(sample_out):
    mol1 = g16.g16_structure(sample_out)
    mol2 = copy.deepcopy(mol1)
    rng = np.random.default_rng(1)
    mol2.xyz = mol1.xyz + rng.normal(scale=0.02, size=mol1.xyz.shape)

    ax = g16.g16_draw_deformation(mol1, mol2, overlay=True, scale=20)
    assert ax is not None


def test_draw_deformation_identical_coords_warns_and_returns_none(sample_out):
    mol1 = g16.g16_structure(sample_out)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ax = g16.g16_draw_deformation(mol1, mol1, scale=20)
    assert ax is None
    assert any("identical coordinates" in str(w.message) for w in caught)


def test_draw_deformation_atom_count_mismatch_raises(sample_out):
    mol1 = g16.g16_structure(sample_out)
    mol_bad = copy.deepcopy(mol1)
    mol_bad.Natoms = mol1.Natoms - 1
    with pytest.raises(ValueError):
        g16.g16_draw_deformation(mol1, mol_bad, scale=20)
