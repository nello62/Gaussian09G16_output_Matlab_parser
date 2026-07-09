import math
import os
import re

import numpy as np

from ._common import Struct, read_lines

_STATE_RE = re.compile(
    r"Excited State\s+(\d+):\s+(\S+)\s+([\d.]+)\s+eV\s+([\d.]+)\s+nm\s+"
    r"f=([\d.]+)(?:\s+<S\*\*2>=([\d.]+))?"
)
_EXC_RE = re.compile(r"^\s*(\d+)\s*->\s*(\d+)\s+(-?[\d.]+)")
_DEEXC_RE = re.compile(r"^\s*(\d+)\s*<-\s*(\d+)\s+(-?[\d.]+)")


def g16_tddft(filename, nstates=math.inf, plot=False, fwhm_ev=0.30, lines=None):
    """Extracts TD-DFT excited states from a Gaussian 16 .out/.log file.

    Parameters
    ----------
    filename : str
    nstates : int, optional — load the first N states only (default: all)
    plot : bool — also generate a UV-Vis plot (default False)
    fwhm_ev : float — Gaussian broadening FWHM in eV (default 0.30)
    lines : list[str], optional

    Returns
    -------
    td : Struct — n, mult, eV, nm, f, S2 (np.ndarray), trans (list of
        (Ncontrib,3) np.ndarray: [MO_from, MO_to, coeff], MO_to<0 for
        de-excitation), Nstates, has_S2, x_nm, eps_cont, FWHM_eV, filename.
    """
    if lines is None:
        lines = read_lines(filename)
    n = len(lines)

    st_n, st_mult, st_eV, st_nm, st_f, st_S2, st_trans = [], [], [], [], [], [], []

    k = 0
    while k < n:
        ln = lines[k]
        m = _STATE_RE.search(ln)
        if m:
            sn = int(m.group(1))
            if sn > nstates:
                k += 1
                continue

            st_n.append(sn)
            st_mult.append(m.group(2))
            st_eV.append(float(m.group(3)))
            st_nm.append(float(m.group(4)))
            st_f.append(float(m.group(5)))
            st_S2.append(float(m.group(6)) if m.group(6) else math.nan)

            trans_block = []
            k += 1
            while k < n:
                ln2 = lines[k]
                m_ex = _EXC_RE.match(ln2)
                m_de = _DEEXC_RE.match(ln2)
                if m_ex:
                    trans_block.append((int(m_ex.group(1)), int(m_ex.group(2)), float(m_ex.group(3))))
                    k += 1
                elif m_de:
                    trans_block.append((int(m_de.group(1)), -int(m_de.group(2)), float(m_de.group(3))))
                    k += 1
                else:
                    break
            st_trans.append(np.array(trans_block) if trans_block else np.zeros((0, 3)))
            continue

        k += 1

    nstates_found = len(st_n)
    if nstates_found == 0:
        raise ValueError(f"g16_tddft: no excited states found in {filename}\n"
                          f"(check that the job uses TD-DFT)")

    st_eV = np.array(st_eV)
    st_f = np.array(st_f)
    st_S2 = np.array(st_S2)

    x_nm = np.arange(100, 1001, 1, dtype=float)
    x_eV_grid = 1239.84193 / x_nm
    sigma_eV = fwhm_ev / (2 * math.sqrt(2 * math.log(2)))

    eps_cont = np.zeros_like(x_nm)
    for s in range(nstates_found):
        if st_f[s] > 0:
            eps_cont += st_f[s] * np.exp(-(x_eV_grid - st_eV[s]) ** 2 / (2 * sigma_eV ** 2))

    td = Struct(
        n=np.array(st_n), mult=st_mult, eV=st_eV, nm=np.array(st_nm), f=st_f, S2=st_S2,
        trans=st_trans, Nstates=nstates_found, has_S2=not np.all(np.isnan(st_S2)),
        x_nm=x_nm, eps_cont=eps_cont, FWHM_eV=fwhm_ev, filename=filename,
    )

    print(f"\n── g16_tddft: {filename} ──")
    print(f"  {nstates_found} excited states read")
    print(f"  {'State':<6}  {'Mult':<14}  {'eV':>8}  {'nm':>8}  {'f':>8}")
    print(f"  {'-'*52}")
    for s in range(min(nstates_found, 20)):
        print(f"  {td.n[s]:<6d}  {td.mult[s]:<14}  {td.eV[s]:8.4f}  {td.nm[s]:8.2f}  {td.f[s]:8.4f}")
    if nstates_found > 20:
        print(f"  ... ({nstates_found} states total)")
    print()

    if plot:
        _plot_tddft(td)

    return td


def _plot_tddft(td):
    import matplotlib.pyplot as plt

    fname = os.path.splitext(os.path.basename(td.filename))[0]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    fig.canvas.manager.set_window_title(fname)

    for s in range(td.Nstates):
        if td.f[s] > 0:
            ax.plot([td.nm[s], td.nm[s]], [0, td.f[s]], color=(0.70, 0.70, 0.70), linewidth=0.9)

    f_max = td.f.max()
    eps_max = td.eps_cont.max()
    scale = f_max / eps_max if eps_max > 0 else 1.0
    ax.plot(td.x_nm, td.eps_cont * scale, color=(0.15, 0.45, 0.80), linewidth=1.5,
            label=f"Gaussian FWHM={td.FWHM_eV:.2f} eV")

    ax.set_xlabel("Wavelength (nm)", fontsize=10)
    ax.set_ylabel("Oscillator strength  f", fontsize=10)
    ax.set_title(fname, fontsize=11)
    ax.legend(loc="upper right", frameon=False)
    ax.set_xlim(200, td.nm.max() * 1.15)

    fig.tight_layout()
    plt.show()
