import G16parser as g16


def test_charges_basic(sample_out):
    ch = g16.g16_charges(sample_out, plot=False)
    assert len(ch.charges) > 0
    assert ch.Natoms == len(ch.charges)
    assert ch.filename == sample_out


def test_charges_dipole_overlay(sample_out):
    ch = g16.g16_charges(sample_out, plot=False, show_dipole=True)
    assert ch.dipole_Debye is not None
    assert ch.dipole_Debye >= 0


def test_route_basic(sample_out):
    route = g16.g16_route(sample_out)
    assert isinstance(route, str)
    assert route.strip().startswith("#")
    # Gaussian wraps the route echo at a fixed column with no regard for
    # word boundaries -- the fixture file has "nosym"/"m cphf=grid=fine"
    # split exactly across a line wrap. Joining continuation lines with an
    # inserted space (the previous behaviour) corrupted this into
    # "nosym m cphf=grid= fine". Assert the real keywords survive intact.
    assert "nosymm" in route
    assert "cphf=grid=fine" in route
    assert "nosym m" not in route
    assert "grid= fine" not in route


def test_restart_route_not_corrupted(tmp_path, sample_out):
    # G16_restart re-wraps the route for the generated .gjf; reading it
    # back must not have merged the last word of one wrapped line into
    # the first word of the next (a regression this specific case caught
    # once already -- the writer's own wrap point looked identical to a
    # Gaussian mid-word wrap from the reader's point of view).
    gjf = g16.g16_restart(sample_out, output=str(tmp_path / "restart.gjf"))
    ginp = g16.g16_read_input(gjf)
    assert "nosymm cphf=grid=fine" in ginp.route


def test_charge_mult_basic(sample_out):
    charge, mult = g16.g16_charge_mult(sample_out)
    assert isinstance(charge, int)
    assert mult >= 1


def test_gaussian_version_basic(sample_out):
    gv = g16.g16_gaussian_version(sample_out)
    assert gv.major is not None
    assert "Gaussian" in gv.full
    assert gv.filename == sample_out


def test_get_bond_length_from_structure(sample_out):
    mol = g16.g16_structure(sample_out)
    bt = g16.g16_get_bond_length(mol, include_h=True)
    bt_noh = g16.g16_get_bond_length(mol, include_h=False)
    assert len(bt) >= len(bt_noh)
