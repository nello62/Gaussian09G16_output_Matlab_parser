import os
import re

import G16parser as g16


def test_split_two_jobs_named_from_wrapped_chk(sample_link1_log, tmp_path):
    files = g16.g16_split_link1_log(sample_link1_log, naming="chk", outdir=str(tmp_path))
    assert len(files) == 2
    name0 = os.path.splitext(os.path.basename(files[0]))[0]
    # The %chk Link0 line of the first job is wrapped across two lines
    # (Gaussian wraps at 80 columns) -- exactly what the continuation-
    # joining logic exists for.
    assert name0 == "NMA_wB97XD_def2TZVP_TD_trans"
    for f in files:
        assert os.path.isfile(f)


def test_split_second_job_readable_and_complete(sample_link1_log, tmp_path):
    files = g16.g16_split_link1_log(sample_link1_log, naming="chk", outdir=str(tmp_path))
    with open(files[1], encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    n_states = len(re.findall(r"^\s*Excited State\s+\d+:", text, re.MULTILINE))
    assert n_states == 40  # 20 singlets + 20 triplets
    assert "Normal termination" in text


def test_split_auto_naming_falls_back_to_chk(sample_link1_log, tmp_path):
    # This fixture's title has no "scan along X: N = ..." pattern, so
    # 'auto' must resolve to 'chk' naming, same result as naming='chk'.
    auto_dir = tmp_path / "auto"
    chk_dir = tmp_path / "chk"
    auto_files = g16.g16_split_link1_log(sample_link1_log, naming="auto", outdir=str(auto_dir))
    chk_files = g16.g16_split_link1_log(sample_link1_log, naming="chk", outdir=str(chk_dir))
    auto_names = [os.path.basename(f) for f in auto_files]
    chk_names = [os.path.basename(f) for f in chk_files]
    assert auto_names == chk_names


def test_split_single_job_file_returns_one_log(sample_out, tmp_path):
    # Every Gaussian log, chained or not, starts with its own "Entering
    # Link 1" line -- a single-job file is simply re-wrapped as one log.
    files = g16.g16_split_link1_log(sample_out, naming="chk", outdir=str(tmp_path))
    assert len(files) == 1


def test_split_field_naming_without_title_raises(sample_link1_log, tmp_path):
    try:
        g16.g16_split_link1_log(sample_link1_log, naming="field", outdir=str(tmp_path))
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_split_non_gaussian_file_raises(tmp_path):
    bogus = tmp_path / "not_a_gaussian_log.log"
    bogus.write_text("this is not a Gaussian output file\n", encoding="utf-8")
    try:
        g16.g16_split_link1_log(str(bogus), outdir=str(tmp_path))
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_split_invalid_naming_raises(sample_link1_log, tmp_path):
    try:
        g16.g16_split_link1_log(sample_link1_log, naming="bogus", outdir=str(tmp_path))
        assert False, "expected ValueError"
    except ValueError:
        pass
