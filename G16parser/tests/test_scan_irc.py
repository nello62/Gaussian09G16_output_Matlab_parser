"""Tests for g16_scan_read / g16_irc_read and the multi-structure warning
of g16_read_all, on small synthetic Gaussian-16-like outputs that contain
only the blocks these parsers read (no real calculation needed)."""

import warnings

import numpy as np
import pytest

import G16parser as g16

SEP = " " + "-" * 69


def _orientation(xyz, label="Input orientation"):
    out = [f"                          {label}:", SEP,
           " Center     Atomic      Atomic             Coordinates (Angstroms)",
           " Number     Number       Type             X           Y           Z", SEP]
    for n, (z, x, y, w) in enumerate(xyz, 1):
        out.append(f" {n:6d} {z:10d} {0:11d} {x:15.6f} {y:11.6f} {w:11.6f}")
    out.append(SEP)
    return out


def _hooh(phi_deg):
    """H-O-O-H with dihedral phi (deg): O1 at origin, O2 on x."""
    p = np.radians(phi_deg)
    return [(1, -0.3, 0.9, 0.0), (8, 0.0, 0.0, 0.0), (8, 1.45, 0.0, 0.0),
            (1, 1.75, 0.9 * np.cos(p), 0.9 * np.sin(p))]


@pytest.fixture
def scan_file(tmp_path):
    L = [" Entering Link 1 = C:\\G16W\\l1.exe PID=      1.", " Gaussian 16, Revision C.01,",
         " The following ModRedundant input section has been read:",
         " D       1       2       3       4 S   2 60.000",
         " B       2       3 F", ""]
    energies = [-151.20, -151.19, -151.17]
    for p, (phi, e) in enumerate(zip([170.0, -130.0, -70.0], energies), 1):
        for step in (1, 2):
            L += _orientation(_hooh(phi + (5.0 if step == 1 and p > 1 else 0.0)))
            L += [f" SCF Done:  E(RwB97XD) =  {e + (0.001 if step == 1 else 0.0):.9f}     A.U. after   9 cycles",
                  f" Step number   {step} out of a maximum of  100 on scan point     {p} out of     3"]
        L += ["    -- Stationary point found."]
    L += [" Normal termination of Gaussian 16 at Mon Nov 02 19:58:00 2026."]
    f = tmp_path / "scan.log"
    f.write_text("\n".join(L) + "\n")
    return str(f)


@pytest.fixture
def irc_file(tmp_path):
    ts = _hooh(90.0)
    L = [" Entering Link 1 = C:\\G16W\\l1.exe PID=      1.", " Gaussian 16, Revision C.01,"]
    L += _orientation(ts) + [" SCF Done:  E(RwB97XD) =  -151.100000000     A.U. after   9 cycles",
                             "    -- Stationary point found.",
                             " Normal termination of Gaussian 16 at Mon Nov 02 21:56:49 2026."]
    L += [" Entering Link 1 = C:\\G16W\\l1.exe PID=      2.",
          " Structure from the checkpoint file:  \"ts.chk\"",
          " Charge =  0 Multiplicity = 1",
          " H,0,-0.3,0.9,0.0", " O,0,0.0,0.0,0.0", " O,0,1.45,0.0,0.0", " H,0,1.75,0.0,0.9",
          " Recover connectivity data from disk.",
          " IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC",
          " Maximum points per path      =  2",
          " Step size                    =   0.100 bohr",
          " Integration scheme           = HPC",
          "   Energy From Chk =   -151.1000000",
          " Point Number:   0          Path Number:   1"]

    def point(n, direction, path, phi, e, rx):
        out = [f" Point Number  {n} in {direction} path direction."]
        out += _orientation(_hooh(phi + 2.0))
        out += [f" SCF Done:  E(RwB97XD) =  {e + 0.0005:.9f}     A.U. after   9 cycles",
                f"   Corrected End Point Energy =  {e:.6f}",
                f" Point Number:  {n}          Path Number:   {path}",
                "                    CURRENT STRUCTURE",
                "                Cartesian Coordinates (Ang):", SEP,
                " Center     Atomic                     Coordinates (Angstroms)",
                " Number     Number                        X           Y           Z", SEP]
        for k, (z, x, y, w) in enumerate(_hooh(phi), 1):
            out.append(f" {k:6d} {z:10d} {x:27.6f} {y:11.6f} {w:11.6f}")
        out += [SEP, f"   NET REACTION COORDINATE UP TO THIS POINT =    {rx:.5f}"]
        return out

    L += point(1, "FORWARD", 1, 100.0, -151.101, 0.1)
    L += point(2, "FORWARD", 1, 110.0, -151.105, 0.2)
    L += [" Beginning calculation of the REVERSE path."]
    L += point(1, "REVERSE", 2, 80.0, -151.102, 0.1)
    L += [" Maximum number of corrector steps exceded.",
          " Error termination via Lnk1e in C:\\G16W\\l123.exe at Tue Nov 03 00:47:37 2026."]
    f = tmp_path / "ts_irc.log"
    f.write_text("\n".join(L) + "\n")
    return str(f)


def test_scan_read(scan_file):
    S = g16.g16_scan_read(scan_file)
    assert [c.label for c in S.coords] == ["D(1,2,3,4)"]
    assert S.coords[0].nsteps == 2 and S.coords[0].stepsize == 60.0
    assert S.frozen == ["B(2,3)"]
    assert S.E.shape == (3,) and S.n_points_expected == 3
    np.testing.assert_allclose(S.E, [-151.20, -151.19, -151.17])
    # dihedral measured on the converged geometries, unwrapped across +-180
    np.testing.assert_allclose(S.value[:, 0], [170.0, 230.0, 290.0], atol=1e-3)
    assert S.i_min == 0 and S.dE_kcal[0] == 0.0
    assert np.isclose(S.dE_kcal[2], 0.03 * 627.5095)
    assert S.converged.all() and S.completed
    assert list(S.opt_steps) == [2, 2, 2]
    assert S.xyz.shape == (3, 4, 3) and S.symbols == ["H", "O", "O", "H"]
    tbl, phi = g16.g16_get_dihedral(S.mol[1], [1, 2, 3, 4])
    assert np.isclose(phi[0], -130.0, atol=1e-3)


def test_scan_read_plot(scan_file):
    S = g16.g16_scan_read(scan_file, plot=True)
    assert S is not None


def test_irc_read(irc_file):
    I = g16.g16_irc_read(irc_file)
    assert I.job == 2
    assert I.direction == ["REVERSE", "TS", "FORWARD", "FORWARD"]
    assert I.i_ts == 1 and I.n_forward == 2 and I.n_reverse == 1
    np.testing.assert_allclose(I.rxcoord, [-0.1, 0.0, 0.1, 0.2])
    np.testing.assert_allclose(I.E, [-151.102, -151.1, -151.101, -151.105])
    assert I.E_source == ["hpc_corrected", "chk", "hpc_corrected", "hpc_corrected"]
    assert I.dE_kcal[I.i_ts] == 0.0 and (I.dE_kcal[[0, 2, 3]] < 0).all()
    assert not I.completed and "corrector" in I.status
    assert I.max_points == 2 and I.scheme == "HPC" and I.step_size == "0.100 bohr"
    # TS geometry from the checkpoint echo; path points from CURRENT STRUCTURE
    _, phi_ts = g16.g16_get_dihedral(I.mol[1], [1, 2, 3, 4])
    _, phi_f2 = g16.g16_get_dihedral(I.mol[3], [1, 2, 3, 4])
    assert np.isclose(phi_ts[0], 90.0, atol=1e-3) and np.isclose(phi_f2[0], 110.0, atol=1e-3)
    assert I.xyz.shape == (4, 4, 3)


def test_irc_read_plot(irc_file):
    assert g16.g16_irc_read(irc_file, plot=True) is not None


def test_wrong_file_type(scan_file, irc_file):
    with pytest.raises(ValueError):
        g16.g16_irc_read(scan_file)
    with pytest.raises(ValueError):
        g16.g16_scan_read(irc_file)


def test_read_all_warns_on_multistructure(irc_file):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        from G16parser.read_all import _warn_multistructure
        from G16parser._common import read_lines
        _warn_multistructure(irc_file, read_lines(irc_file))
    msgs = " ".join(str(x.message) for x in w)
    assert "IRC" in msgs and "--Link1--" in msgs


def test_listed_in_toolbox():
    T = g16.g16_list()
    assert {"g16_scan_read", "g16_irc_read"} <= set(T["Name"])
