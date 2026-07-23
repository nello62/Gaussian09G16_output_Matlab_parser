import os
import shutil

import G16parser as g16


def test_read_all_basic(sample_out):
    T = g16.g16_read_all(sample_out)
    assert T.structure.Natoms > 0
    assert T.energy.SCF < 0
    assert T.dipolar.mu_tot >= 0
    assert T.route.strip().startswith("#")
    assert T.chargemol.mol >= 1


def test_write_report(tmp_path, sample_out):
    T = g16.g16_read_all(sample_out)
    outfile = tmp_path / "report.txt"
    result_path = g16.g16_write_report(T, str(outfile))

    assert result_path == str(outfile)
    assert outfile.exists()
    content = outfile.read_text()
    assert len(content) > 0


def test_restart_roundtrip(tmp_path, sample_out):
    out_gjf = tmp_path / "restarted.gjf"
    gjf_file = g16.g16_restart(sample_out, output=str(out_gjf))

    assert gjf_file == str(out_gjf)
    assert out_gjf.exists()

    mol = g16.g16_structure(sample_out)
    ginp = g16.g16_read_input(gjf_file)
    assert ginp.Natoms == mol.Natoms
    assert set(ginp.symbols) == set(mol.symbols)


def test_batch_read_all(tmp_path, sample_out):
    batch_dir = tmp_path / "batch"
    batch_dir.mkdir()
    shutil.copy(sample_out, batch_dir / os.path.basename(sample_out))

    summary_path = tmp_path / "summary.csv"
    T = g16.g16_batch_read_all(str(batch_dir), save_as=str(summary_path))

    assert len(T) == 1
    assert T.iloc[0]["Status"] == "ok"
    assert T.iloc[0]["Natoms"] > 0
    assert summary_path.exists()


def test_batch_read_all_resilient_to_bad_file(tmp_path, sample_out):
    batch_dir = tmp_path / "batch_mixed"
    batch_dir.mkdir()
    shutil.copy(sample_out, batch_dir / os.path.basename(sample_out))
    (batch_dir / "garbage.log").write_text("this is not a Gaussian output file\n")

    T = g16.g16_batch_read_all(str(batch_dir))

    assert len(T) == 2
    statuses = set(T["Status"])
    assert "ok" in statuses
    assert any(s != "ok" for s in statuses)
