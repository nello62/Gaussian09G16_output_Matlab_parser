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
    being sourced. Excludes files earmarked for a more specific fixture:
    filenames containing "nbo" (case-insensitive, backs sample_nbo_out),
    and any file with a same-stem ".fchk" sibling (backs sample_fchk_pair,
    e.g. tests/fixtures/4-NTP.out/.fchk) -- this fixture backs tests that
    assume a "plain" file (e.g. asserting g16_nbo_bonds raises on a file
    with no NBO section), so it must not resolve to one of those instead.
    """
    candidates = sorted(FIXTURES_DIR.glob("*.out")) + sorted(FIXTURES_DIR.glob("*.log"))
    candidates = [c for c in candidates
                  if "nbo" not in c.name.lower() and not c.with_suffix(".fchk").exists()]
    if not candidates:
        pytest.skip("No Gaussian 16 sample .out/.log file in tests/fixtures/ yet")
    return str(candidates[0])


@pytest.fixture(scope="session")
def sample_nbo_out():
    """Path to a real Gaussian .out/.log file with an NBO analysis
    (pop=nbo), used by g16_nbo_bonds tests. Distinct from sample_out
    since most fixture files (e.g. test.out) have no NBO section.

    Skips (rather than fails) any test that depends on it if no such
    file is present yet.
    """
    candidates = sorted(FIXTURES_DIR.glob("*NBO*")) + sorted(FIXTURES_DIR.glob("*nbo*"))
    if not candidates:
        pytest.skip("No NBO-analysis sample file in tests/fixtures/ yet")
    return str(candidates[0])


@pytest.fixture(scope="session")
def sample_fchk():
    """Path to a real Gaussian formatted checkpoint (.fchk) file dropped
    in tests/fixtures/, used by g16_fchk_read/g16_charges_fchk tests.

    Skips (rather than fails) any test that depends on it if no such
    file is present yet.
    """
    candidates = sorted(FIXTURES_DIR.glob("*.fchk"))
    if not candidates:
        pytest.skip("No .fchk sample file in tests/fixtures/ yet")
    return str(candidates[0])


@pytest.fixture(scope="session")
def sample_fchk_pair():
    """Path pair (fchk, out) for the SAME Gaussian job -- i.e. a
    '<stem>.fchk' with a same-stem '<stem>.out'/'<stem>.log' sibling in
    tests/fixtures/. Used to cross-validate g16_fchk_read's IR/Raman
    intensities (derived from the .fchk dipole/polarisability
    derivatives) against the IR/Raman Gaussian itself prints directly in
    the .out file for g16_nmodes -- these must agree, since both describe
    the same physical quantity for the same calculation.

    Skips (rather than fails) any test that depends on it if no matching
    pair is present yet.
    """
    for fchk in sorted(FIXTURES_DIR.glob("*.fchk")):
        for ext in (".out", ".log"):
            candidate = fchk.with_suffix(ext)
            if candidate.exists():
                return str(fchk), str(candidate)
    pytest.skip("No matching .fchk/.out(.log) pair in tests/fixtures/ yet")


@pytest.fixture(scope="session")
def sample_gjf():
    """Path to the small synthetic Gaussian input file used for
    g16_read_input / g16_restart round-trip tests. Synthetic (hand-written,
    not from a real calculation), so safe to keep committed.
    """
    return str(FIXTURES_DIR / "water.gjf")
