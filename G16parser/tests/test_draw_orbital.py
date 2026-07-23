import G16parser as g16


def test_draw_orbital_basic(sample_out):
    oe = g16.g16_orbital_energies(sample_out)
    ax = g16.g16_draw_orbital(oe)
    assert ax is not None
