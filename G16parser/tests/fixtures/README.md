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

`CH4_NBO.LOG` is a real Gaussian output file (`opt freq=raman
b3lyp/6-311++g(d,p) pop=nbo geom=connectivity polar`, a multi-step
opt+freq+polar job on CH4) used by `test_nbo_bonds.py` via the
`sample_nbo_out` fixture in `conftest.py`. It is the only fixture with an
actual NBO analysis (`pop=nbo`) section, needed to test
`g16_nbo_bonds`; `test.out` deliberately has no NBO data, which is used
to test `g16_nbo_bonds`'s "no NBO analysis found" error path.

`4-NTP.fchk` + `4-NTP.out` are a real, same-molecule pair (opt+freq
B3LYP/6-311G(d,p) on 4-nitrothiophenol, 15 atoms) used by
`test_fchk_read.py` via the `sample_fchk` and `sample_fchk_pair` fixtures
in `conftest.py`. The pair exists specifically to cross-validate
`g16_fchk_read`'s IR/Raman intensities (derived from the `.fchk`
dipole/polarisability derivatives) against the IR/Raman Gaussian itself
prints directly in the `.out` file — this caught a real bug (a data-
scrambling double `reshape` in the MATLAB `G09_fchk_read.m`/
`G16_fchk_read.m`, fixed 2026-08-03; see the toolbox manual's "Notes on
the port" for details) that produced plausible-looking but wrong
IR/Raman values. `sample_out` (used by tests that assume a "plain" file)
explicitly excludes any `.out` with a same-stem `.fchk` sibling, so it
never resolves to `4-NTP.out`.
