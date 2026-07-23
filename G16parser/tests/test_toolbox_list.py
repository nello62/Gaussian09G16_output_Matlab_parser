import G16parser as g16


def test_list_contains_known_functions():
    T = g16.g16_list()
    names = set(T["Name"])

    for expected in (
        "g16_structure", "g16_read_input", "g16_energy", "g16_restart",
        "g16_batch_read_all", "g16_draw_molecule", "g16_write_report",
    ):
        assert expected in names

    assert len(T) >= 20
