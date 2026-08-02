import pytest

import G16parser as g16


def test_nbo_bonds_ch4(sample_nbo_out):
    bt = g16.g16_nbo_bonds(sample_nbo_out)

    assert len(bt) == 4
    assert (bt["BondOrder"] == 1).all()
    assert set(bt["Sym1"]) == {"C"}
    assert set(bt["Sym2"]) == {"H"}


def test_nbo_bonds_first_last_agree(sample_nbo_out):
    bt_last = g16.g16_nbo_bonds(sample_nbo_out, section="last")
    bt_first = g16.g16_nbo_bonds(sample_nbo_out, section="first")

    assert (bt_last["BondOrder"] == bt_first["BondOrder"]).all()
    assert (bt_last["Atom1"] == bt_first["Atom1"]).all()
    assert (bt_last["Atom2"] == bt_first["Atom2"]).all()


def test_nbo_bonds_no_nbo_section_raises(sample_out):
    with pytest.raises(ValueError):
        g16.g16_nbo_bonds(sample_out)


def test_nbo_bonds_double_triple_synthetic(tmp_path):
    synthetic = tmp_path / "synthetic_nbo.log"
    synthetic.write_text(
        " Natural Bond Orbitals (Summary)\n"
        "\n"
        "                                                     Principal Delocalizations\n"
        "           NBO                 Occupancy    Energy   (geminal,vicinal,remote)\n"
        " ===============================================================================\n"
        "     1. BD (   1) C   1 - C   2          1.98000    -0.50000\n"
        "     2. BD (   2) C   1 - C   2          1.97000    -0.49000\n"
        "     3. BD (   1) C   2 - N   3          1.96000    -0.48000\n"
        "     4. BD (   2) C   2 - N   3          1.95000    -0.47000\n"
        "     5. BD (   3) C   2 - N   3          1.94000    -0.46000\n"
        "     6. BD (   1) C   1 - H   4          1.99000    -0.51000\n"
        "     7. BD*(   1) C   1 - C   2          0.02000     0.50000\n"
        "\n"
        "    -------------------------------\n"
        "           Total Lewis\n"
    )

    bt = g16.g16_nbo_bonds(str(synthetic))
    order = {(r.Atom1, r.Atom2): r.BondOrder for _, r in bt.iterrows()}

    assert order[(1, 2)] == 2   # BD* excluded, only the 2 real BD orbitals count
    assert order[(2, 3)] == 3
    assert order[(1, 4)] == 1


def _bond_list_from_nbo(bt):
    bond_list = bt[["Atom1", "Atom2", "BondOrder"]].to_numpy()
    bond_list[:, :2] -= 1   # 1-based (Gaussian) -> 0-based (g16_draw_molecule)
    return bond_list


def test_nbo_bonds_feeds_draw_molecule(sample_nbo_out):
    # Regression: g16_nbo_bonds returns 1-based (Gaussian) atom indices,
    # but g16_draw_molecule's bond_list is 0-based -- passing Atom1/Atom2
    # straight through without the -1 offset previously raised
    # IndexError inside g16_draw_molecule.
    mol = g16.g16_structure(sample_nbo_out)
    bt = g16.g16_nbo_bonds(sample_nbo_out)
    bond_list = _bond_list_from_nbo(bt)

    ax = g16.g16_draw_molecule(mol, bond_list=bond_list)
    assert ax is not None


def test_nbo_bonds_feeds_charges_bond_list(sample_nbo_out):
    bt = g16.g16_nbo_bonds(sample_nbo_out)
    bond_list = _bond_list_from_nbo(bt)

    ch = g16.g16_charges(sample_nbo_out, bond_list=bond_list)
    assert ch.Natoms > 0
