import numpy as np

import G16parser as g16


def test_read_input_water(sample_gjf):
    ginp = g16.g16_read_input(sample_gjf)

    assert ginp.Natoms == 3
    assert ginp.symbols == ["O", "H", "H"]
    assert ginp.xyz.shape == (3, 3)
    assert ginp.charge == 0
    assert ginp.mult == 1
    assert ginp.method.lower() == "b3lyp"
    assert "6-31g(d)" in ginp.basis.lower()
    assert ginp.chk == "water.chk"
    assert ginp.mem == "8GB"
    assert ginp.nproc == "4"
    assert ginp.filename == sample_gjf


def test_read_input_extension_agnostic(tmp_path, sample_gjf):
    content = open(sample_gjf).read()
    for ext in ("com", "in"):
        p = tmp_path / f"water.{ext}"
        p.write_text(content)
        g = g16.g16_read_input(str(p))
        assert g.Natoms == 3
        assert g.charge == 0
        assert g.mult == 1


def test_read_input_compatible_with_draw_molecule(sample_gjf):
    ginp = g16.g16_read_input(sample_gjf)
    ax = g16.g16_draw_molecule(ginp)
    assert ax is not None


def test_read_input_compatible_with_get_bond_length(sample_gjf):
    ginp = g16.g16_read_input(sample_gjf)
    bt = g16.g16_get_bond_length(ginp)
    assert len(bt) == 2  # two O-H bonds
    assert np.allclose(bt["Distance_Ang"].values, 0.9578, atol=1e-3)


def test_read_input_no_geometry_raises(tmp_path):
    p = tmp_path / "bad.gjf"
    p.write_text("#p opt b3lyp/6-31g(d)\n\ntitle\n\n0 1\n\n")
    try:
        g16.g16_read_input(str(p))
        assert False, "expected ValueError"
    except ValueError:
        pass
