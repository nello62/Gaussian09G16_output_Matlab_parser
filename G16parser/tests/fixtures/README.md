# Test fixtures

`water.gjf` is a small hand-written, synthetic Gaussian input file used by
`test_read_input.py` — safe to keep committed (not from a real calculation).

`test.out` is a real Gaussian 16 output file (`opt=calcall freq=raman
field=x-5 CPHF=Rdfreq b3lyp/6-311g(d,p) nosym`) used by `test_structure.py`,
`test_energetics.py`, `test_vibrational.py`, `test_charges_route.py`,
`test_read_all_and_reports.py`, `test_response_properties.py`, and
`test_draw_orbital.py` via the `sample_out` fixture in `conftest.py`. It
covers geometry optimisation (convergence data), vibrational
frequencies/IR/Raman, and static + dynamic hyperpolarisability (Beta) —
it has no TD-DFT excited states, which is fine: `g16_tddft`'s "no excited
states found" case is itself part of what `test_tddft_absent_or_valid`
exercises.

If this file is ever removed, any `.out`/`.log` file dropped into this
folder (any filename) is picked up automatically by the same fixture; the
tests that depend on it are skipped (not failed) when none is present.
