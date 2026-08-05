import G16parser as g16


def test_available_data_basic_shape(sample_out):
    T = g16.g16_available_data(sample_out)
    assert set(T.columns) == {"Function", "Available", "Requires", "Notes"}
    assert "nmodes" in T["Function"].values
    assert "tddft" in T["Function"].values


def test_available_data_matches_actual_extraction(tmp_path):
    # Build a minimal synthetic .out with a real Gaussian route-echo block
    # (the only thing g16_route/g16_available_data reads) for a TD-DFT-only
    # job with no 'freq' keyword, and check the predictions agree with
    # what g16_nmodes/g16_energy actually do on a real no-freq file
    # (the same route as tests/fixtures' sibling scenario in the MATLAB
    # toolbox's Test_TD.out: td=(singlets,nstates=100), no freq).
    content = (
        " ---------------------------------------------------------------\n"
        " # td=(singlets,nstates=100) b3lyp/6-311g(d,p) pop=full\n"
        " ---------------------------------------------------------------\n"
    )
    f = tmp_path / "td_only.out"
    f.write_text(content, encoding="utf-8")

    T = g16.g16_available_data(str(f))
    d = dict(zip(T["Function"], T["Available"]))

    assert d["nmodes"] is False
    assert d["spectra (IR)"] is False
    assert d["spectra (Raman)"] is False
    assert d["energy (thermochemistry: ZPE/H/G/S)"] is False
    assert d["tddft"] is True
    assert d["charges (Mulliken)"] is True
    assert d["hyperpolar (vibrational Beta)"] is False


def test_available_data_freq_raman_job(tmp_path):
    content = (
        " ---------------------------------------------------------------\n"
        " # opt=calcall freq=raman pop=(full,nbo) b3lyp/6-311g(d,p) polar\n"
        " ---------------------------------------------------------------\n"
    )
    f = tmp_path / "freq_raman.out"
    f.write_text(content, encoding="utf-8")

    T = g16.g16_available_data(str(f))
    d = dict(zip(T["Function"], T["Available"]))

    assert d["nmodes"] is True
    assert d["spectra (Raman)"] is True
    assert d["energy (thermochemistry: ZPE/H/G/S)"] is True
    assert d["nbo_bonds"] is True
    assert d["dipole_polar (static Alpha)"] is True
    assert d["dipole_polar (dynamic Alpha)"] is False  # polar present but no cphf=rdfreq
    assert d["convergence"] is True
