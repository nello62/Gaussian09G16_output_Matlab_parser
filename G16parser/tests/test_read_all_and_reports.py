import os
import shutil
import warnings

import G16parser as g16


def test_read_all_basic(sample_out):
    T = g16.g16_read_all(sample_out)
    assert T.structure.Natoms > 0
    assert T.energy.SCF < 0
    assert T.dipolar.mu_tot >= 0
    assert T.route.strip().startswith("#")
    assert T.chargemol.mol >= 1


def test_read_all_no_freq_section_does_not_crash(tmp_path, sample_out):
    # Regression test: a Gaussian job with no 'freq' keyword (e.g. a
    # TD-DFT-only single point job) has no "normal coordinates"/"Harmonic
    # frequencies" section at all, so g16_nmodes/g16_spectra genuinely
    # raise ValueError -- g16_read_all used to let that exception
    # propagate and crash the whole call, even though every other field
    # (charges, energy, structure, dipole, route, charge/mult) would have
    # extracted fine. Simulates a no-freq file by breaking the two
    # section headers in a real fixture, rather than requiring a second
    # large real .out fixture.
    with open(sample_out, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()
    content = content.replace("and normal coordinates:", "and NORMAL_COORDS_REMOVED:")
    content = content.replace("Harmonic frequencies", "HARMONIC_FREQ_REMOVED")

    no_freq_out = tmp_path / "no_freq.out"
    no_freq_out.write_text(content, encoding="utf-8")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        T = g16.g16_read_all(str(no_freq_out))
        messages = [str(w.message) for w in caught]

    assert T.nmodes is None
    assert T.spectra is None
    assert T.charge is not None
    assert T.energy is not None
    assert T.structure.Natoms > 0
    assert T.dipolar is not None
    assert T.route.strip().startswith("#")
    assert any("nmodes omitted" in m for m in messages)
    assert any("spectra omitted" in m for m in messages)

    # g16_write_report must also handle the missing fields gracefully
    report_path = g16.g16_write_report(T, str(tmp_path / "no_freq_report.txt"))
    assert os.path.isfile(report_path)


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
