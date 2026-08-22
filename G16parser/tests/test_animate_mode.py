import json
import shutil
import subprocess

import pytest

import G16parser as g16

pytestmark = pytest.mark.skipif(shutil.which("ffmpeg") is None, reason="ffmpeg not on PATH")


def _video_dims(path):
    result = subprocess.run(
        ["ffprobe", "-v", "error", "-select_streams", "v:0",
         "-show_entries", "stream=width,height", "-of", "json", str(path)],
        capture_output=True, text=True, check=True,
    )
    stream = json.loads(result.stdout)["streams"][0]
    return stream["width"], stream["height"]


def test_animate_mode_crop_is_smaller_than_uncropped(sample_nbo_out, tmp_path):
    mol = g16.g16_structure(sample_nbo_out)
    nm = g16.g16_nmodes(sample_nbo_out)

    out_full = tmp_path / "full.mp4"
    out_crop = tmp_path / "crop.mp4"

    g16.g16_animate_mode(mol, nm, 1, filename=str(out_full), crop=False,
                          frames_per_cycle=6, n_cycles=1, fps=10)
    g16.g16_animate_mode(mol, nm, 1, filename=str(out_crop), crop=True,
                          frames_per_cycle=6, n_cycles=1, fps=10)

    assert out_full.exists() and out_crop.exists()

    w_full, h_full = _video_dims(out_full)
    w_crop, h_crop = _video_dims(out_crop)

    assert w_crop < w_full and h_crop < h_full
    assert w_crop % 2 == 0 and h_crop % 2 == 0  # required by yuv420p


def test_animate_mode_show_title_enlarges_crop_bbox(sample_nbo_out, tmp_path):
    mol = g16.g16_structure(sample_nbo_out)
    nm = g16.g16_nmodes(sample_nbo_out)

    out_notitle = tmp_path / "notitle.mp4"
    out_title = tmp_path / "title.mp4"

    g16.g16_animate_mode(mol, nm, 1, filename=str(out_notitle), crop=True,
                          show_title=False, frames_per_cycle=6, n_cycles=1, fps=10)
    g16.g16_animate_mode(mol, nm, 1, filename=str(out_title), crop=True,
                          show_title=True, frames_per_cycle=6, n_cycles=1, fps=10)

    w0, h0 = _video_dims(out_notitle)
    w1, h1 = _video_dims(out_title)

    # The title text adds ink above the molecule, so its crop bbox should
    # be at least as tall (never smaller) than the no-title one.
    assert h1 >= h0
