import numpy as np

import G16parser as g16


def test_structure_basic(sample_out):
    mol = g16.g16_structure(sample_out)

    assert mol.Natoms > 0
    assert len(mol.symbols) == mol.Natoms
    assert mol.xyz.shape == (mol.Natoms, 3)
    assert mol.Z.shape == (mol.Natoms,)
    assert mol.n_steps >= 1
    assert 1 <= mol.step <= mol.n_steps
    assert mol.filename == sample_out
    assert not np.isnan(mol.xyz).any()


def test_structure_first_vs_last_step(sample_out):
    mol_last = g16.g16_structure(sample_out, step="last")
    if mol_last.n_steps < 2:
        return  # single-step file: nothing to compare
    mol_first = g16.g16_structure(sample_out, step="first")
    assert mol_first.step == 1
    assert mol_last.step == mol_last.n_steps
    # different optimisation steps should not be bit-identical geometries
    assert not np.allclose(mol_first.xyz, mol_last.xyz)


def test_structure_compatible_with_draw_molecule(sample_out):
    mol = g16.g16_structure(sample_out)
    ax = g16.g16_draw_molecule(mol)
    assert ax is not None


def test_structure_compatible_with_get_bond_length(sample_out):
    mol = g16.g16_structure(sample_out)
    bt = g16.g16_get_bond_length(mol)
    assert len(bt) > 0
    assert (bt["Distance_Ang"].values > 0).all()
