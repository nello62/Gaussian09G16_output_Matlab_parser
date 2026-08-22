import tkinter as tk
from tkinter import messagebox, ttk

import numpy as np

from G16parser import g16_draw_mode, g16_nmodes, g16_structure


def _running_ipython_tk_loop():
    """True if IPython (e.g. Spyder's console) is already pumping a Tk
    event loop, in which case we must NOT call root.mainloop() ourselves
    (it would block the console instead of returning to the prompt).
    Same helper as G16parser's own g16_mode_viewer, duplicated locally
    rather than importing a leading-underscore (private) symbol from a
    foreign package's internal module.
    """
    try:
        ip = get_ipython()  # noqa: F821
    except NameError:
        return False
    return getattr(ip, "active_eventloop", None) == "tk"


def g16_raman_browser(filename, property="raman", scale=1.5, atom_scale=0.35, show_labels=False):
    """Interactive click-to-select stick spectrum browser for vibrational
    modes: plots intensity vs. frequency as a stem ("stick") plot,
    clicking a stick selects the nearest mode, then a button renders its
    3D displacement structure.

    Python port of G_raman_browser.m. This is a standalone G_Utility
    companion to G16parser's own g16_mode_viewer: that tool selects a
    mode from a drop-down list; this one selects it by clicking directly
    on its stick in the spectrum.

    G16-only (see this package's docstring for why -- G16parser itself
    has no Gaussian 09 support, unlike the MATLAB original).

    Parameters
    ----------
    filename : str -- a Gaussian 16 output file (.log/.out).
    property : 'raman' (default) | 'ir' -- which stick intensity to
        plot. 'raman' requires the file to actually contain Raman
        activities (a Freq=Raman job); a clear error suggests 'ir'
        instead otherwise.
    scale, atom_scale, show_labels : forwarded to g16_draw_mode.

    Click anywhere in the spectrum axes to select the nearest mode by
    frequency (highlighted in red); the info line above the plot updates
    immediately. Press "Draw mode" to open (or replace) a 3D structure
    figure showing that mode's displacement arrows.

    If the file has CPHF=RdFreq (pre-resonance Raman) data --
    nm.has_RamanFr, once g16_nmodes exposes it -- the info line reports
    both the static and each dynamic Raman activity, mirroring
    G_raman_browser.m's current behaviour; today (G16parser's g16_nmodes
    does not yet have .RamanFr/.IncidentLight/.has_RamanFr) it falls back
    to the plain single-value form automatically and needs no changes
    once that field is ported.

    Example
    -------
    g16_raman_browser('4-NTP.out')
    g16_raman_browser('4-NTP.out', property='ir', scale=2)
    """
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure

    prop = property.lower()
    if prop not in ("raman", "ir"):
        raise ValueError("g16_raman_browser: property must be 'raman' or 'ir'.")

    print(f"g16_raman_browser: reading structure and normal modes from {filename} ...")
    mol = g16_structure(filename)
    nm = g16_nmodes(filename)
    print(f"  {nm.Natoms} atoms, {nm.Nmodes} vibrational modes.")

    if nm.Nmodes < 1:
        raise ValueError(f"g16_raman_browser: no vibrational modes found in {filename}.")
    if prop == "raman" and not nm.has_Raman:
        raise ValueError(
            f"g16_raman_browser: {filename} has no Raman activities (not a Freq=Raman job) -- "
            "pass property='ir' to browse the IR spectrum instead."
        )

    vals = nm.Raman if prop == "raman" else nm.IR
    ylab = "Raman activity (Å^4/AMU)" if prop == "raman" else "IR intensity (KM/Mole)"
    freqs = np.asarray(nm.freq)
    has_ramanfr = prop == "raman" and getattr(nm, "has_RamanFr", False)

    state = {"selected_idx": None}   # 1-based mode index (matches g16_draw_mode)

    # ---- build the browser window -----------------------------------------
    # Reuse an already-existing Tk root/interpreter if one is active, same
    # rationale as g16_mode_viewer: a second independent Tcl/Tk interpreter
    # in the same process is a known crash source on macOS.
    _existing_root = tk._default_root
    standalone_root = _existing_root is None
    root = tk.Tk() if standalone_root else tk.Toplevel(_existing_root)
    root.title(f"g16_raman_browser - {filename}")
    root.geometry("820x600+100+80")

    frm = ttk.Frame(root, padding=10)
    frm.pack(fill="both", expand=True)

    info_var = tk.StringVar(value=f"Click a stick to select a mode ({nm.Nmodes} modes, {prop}).")
    ttk.Label(frm, textvariable=info_var, font=("", 11, "bold")).pack(anchor="w")

    fig = Figure(figsize=(8, 4.6))
    ax = fig.add_subplot(111)
    markerline, stemlines, baseline = ax.stem(freqs, vals, markerfmt="o", basefmt=" ")
    markerline.set_markersize(4)
    markerline.set_color((0.20, 0.40, 0.80))
    stemlines.set_color((0.20, 0.40, 0.80))
    stemlines.set_linewidth(1.1)
    (highlight,) = ax.plot([], [], "o", markersize=10, markerfacecolor=(0.85, 0.10, 0.10),
                           markeredgecolor=(0.85, 0.10, 0.10))
    ax.set_xlabel("Frequency (cm$^{-1}$)")
    ax.set_ylabel(ylab)
    ax.grid(True)

    canvas = FigureCanvasTkAgg(fig, master=frm)
    canvas.get_tk_widget().pack(fill="both", expand=True, pady=(5, 10))

    btn_row = ttk.Frame(frm)
    btn_row.pack(fill="x")
    draw_btn = ttk.Button(btn_row, text="Draw mode")
    draw_btn.pack(side="left")
    ttk.Label(btn_row, text=f"File: {filename}", foreground="#737373").pack(side="left", padx=(15, 0))

    # ---- callbacks ----------------------------------------------------------
    def on_click(event):
        if event.inaxes is not ax or event.xdata is None:
            return
        k = int(np.argmin(np.abs(freqs - event.xdata)))   # 0-based
        state["selected_idx"] = k + 1                     # 1-based, matches g16_draw_mode
        highlight.set_data([freqs[k]], [vals[k]])
        canvas.draw_idle()

        sym_label = f"  ({nm.symmetry[k]})" if nm.symmetry and nm.symmetry[k] else ""
        if has_ramanfr:
            fr_text = ""
            for c in range(1, nm.RamanFr.shape[1]):
                wl_cm1 = nm.IncidentLight[c]
                if wl_cm1 > 0:
                    wl_nm = 1e7 / wl_cm1
                    fr_text += f", {wl_nm:.0f}nm = {nm.RamanFr[k, c]:.2f}"
            info_var.set(f"Selected: mode {k+1}, {freqs[k]:.1f} cm-1{sym_label}, "
                          f"Raman static = {nm.RamanFr[k, 0]:.2f}{fr_text}")
        else:
            info_var.set(f"Selected: mode {k+1}, {freqs[k]:.1f} cm-1{sym_label}, {prop} = {vals[k]:.2f}")

    def draw_selected_mode():
        k = state["selected_idx"]
        if k is None:
            messagebox.showinfo("No mode selected", "Click a stick in the spectrum first to select a mode.",
                                 parent=root)
            return
        try:
            ax_mode = g16_draw_mode(mol, nm, k, scale=scale, atom_scale=atom_scale, show_labels=show_labels)
        except Exception as e:
            messagebox.showerror("g16_draw_mode error", str(e), parent=root)
            return
        mode_fig = ax_mode.figure
        mode_fig.canvas.manager.set_window_title(f"{filename} - Mode {k} ({freqs[k-1]:.1f} cm-1)")
        mode_fig.canvas.manager.show()

    canvas.mpl_connect("button_press_event", on_click)
    draw_btn.config(command=draw_selected_mode)

    def on_close():
        plt.close("all")
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_close)

    if not _running_ipython_tk_loop():
        root.mainloop()
