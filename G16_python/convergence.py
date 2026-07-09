import math
import re

import numpy as np

from ._common import Struct, read_lines, fortran_to_float

_HEADER = "Item               Value     Threshold  Converged?"
_LINE_RE = re.compile(r"([\d.]+(?:D[+-]?\d+)?)\s+([\d.]+(?:D[+-]?\d+)?)\s+(YES|NO)")


def _parse_conv_line(ln):
    m = _LINE_RE.search(ln)
    if not m:
        return None
    v1 = fortran_to_float(m.group(1))
    v2 = fortran_to_float(m.group(2))
    yes = m.group(3) == "YES"
    return v1, v2, yes


def g16_convergence(filename, plot=False, lines=None):
    """Extracts geometry-optimisation convergence criteria from a Gaussian
    16 .out/.log file.

    Parameters
    ----------
    filename : str
    plot : bool — show a 2x2 semilogy convergence plot (default False)
    lines : list[str], optional — pre-read file lines (see g16_read_all)

    Returns
    -------
    cv : Struct with fields MaxForce, RMSForce, MaxDisp, RMSDisp (np.ndarray),
        thr_MaxForce/thr_RMSForce/thr_MaxDisp/thr_RMSDisp, converged,
        conv_step, Nsteps, filename.
    """
    if lines is None:
        lines = read_lines(filename)
    n = len(lines)

    max_force, rms_force, max_disp, rms_disp = [], [], [], []
    thr_mxf = thr_rmsf = thr_mxd = thr_rmsd = math.nan
    converged = False
    conv_step = math.nan
    step_count = 0

    k = 0
    while k < n:
        if _HEADER in lines[k]:
            if k + 4 >= n:
                break
            rows = [_parse_conv_line(lines[k + i]) for i in (1, 2, 3, 4)]
            if all(r is not None for r in rows):
                step_count += 1
                max_force.append(rows[0][0])
                rms_force.append(rows[1][0])
                max_disp.append(rows[2][0])
                rms_disp.append(rows[3][0])

                if math.isnan(thr_mxf):
                    thr_mxf, thr_rmsf, thr_mxd, thr_rmsd = (
                        rows[0][1], rows[1][1], rows[2][1], rows[3][1]
                    )

                all_yes = rows[0][2] and rows[1][2] and rows[2][2] and rows[3][2]
                if all_yes and not converged:
                    converged = True
                    conv_step = step_count
            k += 5
            continue
        k += 1

    cv = Struct(
        MaxForce=np.array(max_force), RMSForce=np.array(rms_force),
        MaxDisp=np.array(max_disp), RMSDisp=np.array(rms_disp),
        thr_MaxForce=thr_mxf, thr_RMSForce=thr_rmsf,
        thr_MaxDisp=thr_mxd, thr_RMSDisp=thr_rmsd,
        converged=converged, conv_step=conv_step,
        Nsteps=step_count, filename=filename,
    )

    print(f"\ng16_convergence: {filename}")
    print(f"  Steps read  : {step_count}")
    if converged:
        print(f"  Converged   : YES  (step {conv_step})")
    else:
        print("  Converged   : NO")
    if step_count > 0:
        print(f"  Last step   : MaxF={max_force[-1]:.2e} (thr {thr_mxf:.2e})  "
              f"RMSF={rms_force[-1]:.2e} (thr {thr_rmsf:.2e})")
        print(f"                MaxD={max_disp[-1]:.2e} (thr {thr_mxd:.2e})  "
              f"RMSD={rms_disp[-1]:.2e} (thr {thr_rmsd:.2e})")
    print()

    if plot and step_count > 0:
        _plot_convergence(cv)

    return cv


def _plot_convergence(cv):
    import os
    import matplotlib.pyplot as plt

    fname = os.path.splitext(os.path.basename(cv.filename))[0]
    steps = np.arange(1, cv.Nsteps + 1)

    colors = ["#2673CC", "#D9331F", "#1AA64C", "#E68C00"]
    labels = ["Max Force", "RMS Force", "Max Displacement", "RMS Displacement"]
    data_all = [cv.MaxForce, cv.RMSForce, cv.MaxDisp, cv.RMSDisp]
    thrs = [cv.thr_MaxForce, cv.thr_RMSForce, cv.thr_MaxDisp, cv.thr_RMSDisp]

    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    fig.canvas.manager.set_window_title(fname)

    for i, ax in enumerate(axes.flat):
        ax.semilogy(steps, data_all[i], color=colors[i], linewidth=1.5,
                    marker="o", markersize=4, markerfacecolor=colors[i])
        ax.axhline(thrs[i], linestyle="--", color="black", linewidth=1.0)
        if cv.converged and cv.conv_step <= cv.Nsteps:
            ax.semilogy(cv.conv_step, data_all[i][cv.conv_step - 1], marker="*",
                        markersize=12, markerfacecolor="#33CC33", markeredgecolor="black",
                        linestyle="none")
        ax.grid(True, which="both", axis="y")
        ax.set_xlabel("Opt step", fontsize=9)
        ax.set_ylabel(labels[i], fontsize=9)
        ax.set_title(labels[i], fontsize=10)
        ax.set_xlim(1, max(steps))

    fig.suptitle(fname.replace("_", r"\_"), fontsize=11)
    fig.tight_layout()
    plt.show()
