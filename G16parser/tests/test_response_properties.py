import math

import G16parser as g16


def test_hyperpolar_does_not_crash(sample_out):
    hp = g16.g16_hyperpolar(sample_out)
    assert hp.filename == sample_out
    # data may or may not be present in a given sample file; either way it
    # must not raise, and beta_vec must be a float (possibly NaN)
    assert isinstance(hp.beta0.beta_vec, float)


def test_tddft_absent_or_valid(sample_out):
    try:
        td = g16.g16_tddft(sample_out)
    except ValueError:
        return  # sample file has no TD-DFT excited states: acceptable
    assert len(td.eV) > 0
    assert (td.eV > 0).all()


def test_hyperpolar_invalid_units_raises(sample_out):
    try:
        g16.g16_hyperpolar(sample_out, units="bogus")
        assert False, "expected ValueError"
    except ValueError:
        pass
