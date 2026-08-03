import os

import numpy as np

from ._common import Struct


def g16_spectra_nm(nm, FWHM=10, xmin=0, xmax=4000, dx=1, normalize=False, plot=False):
    """Generates Lorentzian-broadened IR/Raman spectra from an
    already-parsed nm (normal modes) struct, without reading or
    re-parsing a Gaussian output file.

    Accepts the nm Struct returned by g16_nmodes, or the .nm sub-struct
    returned by g16_fchk_read — both share the same field layout (freq,
    IR, Raman, has_Raman, Nmodes, filename), so this single function
    works with either regardless of data source (.out/.log vs .fchk).
    This fills a gap left by g16_fchk_read: its output Struct's .nm field
    has frequencies and IR/Raman intensities, but — unlike g16_spectra —
    does not itself build the broadened continuum spectrum.

    Same broadening algorithm and parameters as g16_spectra (peak-
    normalised Lorentzian convolution), just applied to an already-parsed
    nm struct instead of reading a file.

    Parameters
    ----------
    nm : Struct — must have freq, IR, Nmodes; Raman/has_Raman/filename
        optional.
    FWHM, xmin, xmax, dx, normalize, plot : see g16_spectra.

    Returns
    -------
    sp : Struct — freq, IR, Raman (np.ndarray, empty if absent), Nmodes,
        has_Raman, x, IR_cont, Raman_cont, FWHM, filename.

    Example
    -------
        data = g16.g16_fchk_read('molecule.fchk')
        sp = g16.g16_spectra_nm(data.nm, FWHM=15, plot=True)

        nm = g16.g16_nmodes('molecule.out')
        sp = g16.g16_spectra_nm(nm, xmin=400, normalize=True)
    """
    for field in ("freq", "IR", "Nmodes"):
        if not hasattr(nm, field):
            raise ValueError(
                f'g16_spectra_nm: nm is missing field "{field}". Use the nm struct '
                "returned by g16_nmodes, or the nm field from g16_fchk_read."
            )

    freqs = np.asarray(nm.freq, dtype=float).ravel()
    IRs = np.asarray(nm.IR, dtype=float).ravel()
    nmodes = nm.Nmodes

    if nmodes == 0 or freqs.size == 0:
        raise ValueError("g16_spectra_nm: nm contains no vibrational modes (Nmodes = 0).")

    raman = getattr(nm, "Raman", None)
    has_raman = raman is not None and len(raman) == nmodes
    Ramans = np.asarray(raman, dtype=float).ravel() if has_raman else np.array([])

    filename = getattr(nm, "filename", "") or ""

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
        _plot_spectra_nm(sp)

    return sp


def _plot_spectra_nm(sp):
    import matplotlib.pyplot as plt

    fname = os.path.splitext(os.path.basename(sp.filename))[0] if sp.filename else "nm data"
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
