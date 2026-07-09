import os
import re

import numpy as np

from ._common import Struct, read_lines

_RAMAN_SEC_RE = re.compile(r"Harmonic frequencies.*Raman scattering")
_IR_SEC_RE = re.compile(r"Harmonic frequencies.*IR intensities.*KM")
_END_RE = re.compile(r"^\s*(-{10,}|Thermochemistry|Zero-point|Normal termination)")
_FREQ_RE = re.compile(r"^\s*Frequencies\s*--")
_IR_RE = re.compile(r"^\s*IR Inten\s*--")
_RAMAN_RE = re.compile(r"^\s*Raman Activ\s*--")


def _parse_after_dashdash(ln):
    idx = ln.find("--")
    return ln[idx + 2:] if idx != -1 else ln


def _floats(s):
    return [float(x) for x in s.split()]


def g16_spectra(filename, FWHM=10, xmin=0, xmax=4000, dx=1, normalize=False,
                plot=False, section="last", lines=None):
    """Extracts IR and Raman spectra from a Gaussian 16 .out/.log file and
    generates Lorentzian-broadened continuous spectra.

    Parameters mirror G16_spectra.m — see its docstring for details.

    Returns
    -------
    sp : Struct — freq, IR, Raman (np.ndarray, empty if absent), Nmodes,
        has_Raman, x, IR_cont, Raman_cont, FWHM, filename.
    """
    if lines is None:
        lines = read_lines(filename)

    idx_raman_sec = [i for i, ln in enumerate(lines) if _RAMAN_SEC_RE.search(ln)]
    idx_ir_sec = [i for i, ln in enumerate(lines) if _IR_SEC_RE.search(ln)]
    all_sec = sorted(set(idx_raman_sec) | set(idx_ir_sec))

    if not all_sec:
        raise ValueError(f'g16_spectra: no "Harmonic frequencies" section found in {filename}')

    sec_req = section.lower()
    if sec_req == "last":
        sec_start = all_sec[-1]
    elif sec_req == "first":
        sec_start = all_sec[0]
    else:
        raise ValueError("g16_spectra: 'section' must be 'last' or 'first'.")

    freqs, IRs, Ramans = [], [], []
    for k in range(sec_start, len(lines)):
        ln = lines[k]
        if k > sec_start and _END_RE.match(ln):
            break
        if _FREQ_RE.match(ln):
            freqs.extend(_floats(_parse_after_dashdash(ln)))
            continue
        if _IR_RE.match(ln):
            IRs.extend(_floats(_parse_after_dashdash(ln)))
            continue
        if _RAMAN_RE.match(ln):
            Ramans.extend(_floats(_parse_after_dashdash(ln)))
            continue

    nmodes = len(freqs)
    if nmodes == 0:
        raise ValueError("g16_spectra: no frequencies read from the selected section.")

    freqs = np.array(freqs, dtype=float)
    IRs = np.array(IRs[:nmodes], dtype=float)
    Ramans = np.array(Ramans[:nmodes], dtype=float)
    if IRs.size < nmodes:
        IRs = np.concatenate([IRs, np.zeros(nmodes - IRs.size)])

    has_raman = Ramans.size == nmodes and nmodes > 0

    x = np.arange(xmin, xmax + dx / 2, dx)
    gamma = FWHM / 2
    IR_cont = np.zeros_like(x)
    Raman_cont = np.zeros_like(x)
    for m in range(nmodes):
        L = (gamma ** 2) / ((x - freqs[m]) ** 2 + gamma ** 2)
        IR_cont += IRs[m] * L
        if has_raman:
            Raman_cont += Ramans[m] * L

    if normalize:
        if IR_cont.max() > 0:
            IR_cont = IR_cont / IR_cont.max()
        if has_raman and Raman_cont.max() > 0:
            Raman_cont = Raman_cont / Raman_cont.max()

    sp = Struct(
        freq=freqs, IR=IRs, Raman=(Ramans if has_raman else np.array([])),
        has_Raman=has_raman, Nmodes=nmodes, x=x,
        IR_cont=IR_cont, Raman_cont=(Raman_cont if has_raman else np.array([])),
        FWHM=FWHM, filename=filename,
    )

    if plot:
        _plot_spectra(sp)

    return sp


def _plot_spectra(sp):
    import matplotlib.pyplot as plt

    fname = os.path.splitext(os.path.basename(sp.filename))[0]
    nrows = 2 if sp.has_Raman else 1

    fig, axes = plt.subplots(nrows, 1, figsize=(7, 4 * nrows))
    fig.canvas.manager.set_window_title(fname)
    axes = np.atleast_1d(axes)

    row = 0
    if sp.has_Raman:
        ax1 = axes[row]; row += 1
        for m in range(sp.Nmodes):
            if sp.Raman[m] > 0:
                ax1.plot([sp.freq[m], sp.freq[m]], [0, sp.Raman[m]],
                         color=(0.75, 0.75, 0.75), linewidth=0.8)
        ax1.plot(sp.x, sp.Raman_cont, color=(0.15, 0.45, 0.80), linewidth=1.5,
                 label=f"Raman (FWHM = {sp.FWHM:g} cm-1)")
        ax1.invert_xaxis()
        ax1.set_xlabel("Wavenumber (cm-1)", fontsize=10)
        ax1.set_ylabel("Raman activity (A^4 AMU^-1)", fontsize=10)
        ax1.set_title(f"Raman - {fname}", fontsize=11)
        ax1.legend(loc="upper right", frameon=False)
        ax1.set_xlim(sp.x[0], sp.x[-1])

    ax2 = axes[row]
    for m in range(sp.Nmodes):
        if sp.IR[m] > 0:
            ax2.plot([sp.freq[m], sp.freq[m]], [0, sp.IR[m]],
                     color=(0.75, 0.75, 0.75), linewidth=0.8)
    ax2.plot(sp.x, sp.IR_cont, color=(0.85, 0.20, 0.15), linewidth=1.5,
             label=f"IR (FWHM = {sp.FWHM:g} cm-1)")
    ax2.invert_xaxis()
    ax2.set_xlabel("Wavenumber (cm-1)", fontsize=10)
    ax2.set_ylabel("IR intensity (KM mol^-1)", fontsize=10)
    ax2.set_title(f"IR - {fname}", fontsize=11)
    ax2.legend(loc="upper right", frameon=False)
    ax2.set_xlim(sp.x[0], sp.x[-1])

    fig.tight_layout()
    plt.show()
