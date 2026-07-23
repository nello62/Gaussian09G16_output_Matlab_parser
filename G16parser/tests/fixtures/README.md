# Test fixtures

`water.gjf` is a small hand-written, synthetic Gaussian input file used by
`test_read_input.py` — safe to keep committed (not from a real calculation).

Most of the test suite (`test_structure.py`, `test_energetics.py`,
`test_vibrational.py`, `test_charges_route.py`,
`test_read_all_and_reports.py`, `test_response_properties.py`,
`test_draw_orbital.py`) needs one real Gaussian 16 `.out`/`.log` file to
run against. Drop any such file in this folder (any filename, `.out` or
`.log` extension) and those tests will pick it up automatically via the
`sample_out` fixture in `conftest.py`. Until a file is present, those
tests are skipped (not failed) so the suite stays green.

For full coverage, an ideal file is an `opt freq` (or `opt freq=raman`)
job so the vibrational/convergence tests exercise real data instead of
just their "no data present" fallback paths.
