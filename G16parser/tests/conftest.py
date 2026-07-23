import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (import before G16parser: see below)

# G16parser (see src/G16parser/_common.py) only force-switches the backend
# to TkAgg if matplotlib.pyplot has not been imported yet by the time it is
# imported. Importing pyplot here first (with Agg already selected) means
# G16parser sees pyplot already loaded and leaves the backend alone instead
# of forcing Tk, which would break headless/CI test runs.
warnings.filterwarnings("ignore", message="matplotlib.pyplot was imported before G16parser")

import pytest  # noqa: E402

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="session")
def sample_out():
    """Path to a real Gaussian 16 .out/.log file dropped in tests/fixtures/.

    Skips (rather than fails) any test that depends on it if no such file
    is present yet, so the suite stays runnable while a fixture file is
    being sourced.
    """
    candidates = sorted(FIXTURES_DIR.glob("*.out")) + sorted(FIXTURES_DIR.glob("*.log"))
    if not candidates:
        pytest.skip("No Gaussian 16 sample .out/.log file in tests/fixtures/ yet")
    return str(candidates[0])


@pytest.fixture(scope="session")
def sample_gjf():
    """Path to the small synthetic Gaussian input file used for
    g16_read_input / g16_restart round-trip tests. Synthetic (hand-written,
    not from a real calculation), so safe to keep committed.
    """
    return str(FIXTURES_DIR / "water.gjf")
