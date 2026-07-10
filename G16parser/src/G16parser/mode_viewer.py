import tkinter as tk
from tkinter import colorchooser, filedialog, messagebox, ttk

import numpy as np

from .draw_mode import g16_draw_mode
from .nmodes import g16_nmodes
from .structure import g16_structure

_COLOR_PRESETS = {
    "Orange (default)": (1.0, 0.4, 0.1),
    "Red": (0.85, 0.0, 0.0),
    "Blue": (0.0, 0.35, 0.85),
    "Green": (0.10, 0.60, 0.20),
    "Black": (0.0, 0.0, 0.0),
    "Magenta": (0.75, 0.0, 0.75),
}
_COLOR_NAMES = list(_COLOR_PRESETS) + ["Custom..."]
_FMT_TO_EXT = {"PDF (vector)": ".pdf", "EPS (vector)": ".eps", "JPEG (raster)": ".jpg"}


def _to_hex(rgb):
    return "#%02x%02x%02x" % tuple(int(round(c * 255)) for c in rgb)


def _running_ipython_tk_loop():
    """True if IPython (e.g. Spyder's console) is already pumping a Tk
    event loop, in which case we must NOT call root.mainloop() ourselves
    (it would block the console instead of returning to the prompt)."""
    try:
        ip = get_ipython()  # noqa: F821
    except NameError:
        return False
    return getattr(ip, "active_eventloop", None) == "tk"


def g16_mode_viewer(filename, scale=1.5, atom_scale=0.35, bond_tol=1.30,
                     show_labels=False, flip_sign=False, arrow_color=None,
                     font_size=10, **extra_opts):
    """Interactive Tkinter selector window for Gaussian vibrational modes.

    Python port of G16_modeViewer.m: reads the structure and all
    vibrational normal modes via g16_structure/g16_nmodes, then opens a
    selector window listing every mode (index, frequency, symmetry).
    Choosing an entry draws that mode with g16_draw_mode in a new
    matplotlib window; previous mode windows are left open so different
    modes can be compared side by side. An "Order by" box re-sorts the
    list by mode number, IR intensity, or Raman intensity (if present).
    Every g16_draw_mode option (Scale, ArrowColor, AtomScale, BondTol,
    ShowLabels, FlipSign) has a matching control that redraws the current
    mode in place when changed, plus a title box and a Save-figure button
    (PDF/EPS/JPEG).

    Simplification vs. the MATLAB original: "target figure" resolution
    always uses the most recently drawn mode window (Tkinter/matplotlib
    has no simple equivalent of MATLAB's groot().CurrentFigure to track
    whichever window the user last clicked on across separate widgets).

    Known issue: running this inside Spyder's IPython console (or any
    ipykernel with the Graphics backend set to Tkinter/Automatic) can
    segfault the whole kernel — a crash inside ipykernel's own periodic
    Tk event-loop pump (eventloops.loop_tk), not in this code; confirmed
    to run correctly as a plain script/terminal session. Prefer:
        python3 -c "import G16parser as g16; g16.g16_mode_viewer('file.out')"
    over calling it from Spyder's console until upstream ipykernel/Tk
    support for a full custom Tkinter app (beyond matplotlib's own figure
    windows) is more solid.

    Parameters
    ----------
    filename : str
    scale, atom_scale, bond_tol, show_labels, flip_sign, arrow_color :
        initial values for the matching g16_draw_mode option controls
        (arrow_color=None selects the "Orange (default)" preset; pass an
        (r, g, b) tuple to start on "Custom...").
    font_size : int — initial atom-label font size (post-processing only,
        applied directly to the rendered label Text objects; only visible
        when ShowLabels is on).
    **extra_opts : forwarded as-is to every g16_draw_mode call.
    """
    print(f"g16_mode_viewer: reading structure and normal modes from {filename} ...")
    mol = g16_structure(filename)
    nm = g16_nmodes(filename)
    print(f"  {nm.Natoms} atoms, {nm.Nmodes} vibrational modes.")

    if nm.Nmodes < 1:
        raise ValueError(f"g16_mode_viewer: no vibrational modes found in {filename}.")

    state = {
        "current_fig": None,
        "open_figs": [],
        "custom_color": tuple(arrow_color) if arrow_color is not None else (1.0, 0.4, 0.1),
    }

    def build_mode_items(order_name):
        if order_name == "IR intensity":
            idx = list(np.argsort(-nm.IR))
        elif order_name == "Raman intensity":
            idx = list(np.argsort(-nm.Raman))
        else:
            idx = list(range(nm.Nmodes))

        items = []
        for i in idx:
            sym_label = f"  ({nm.symmetry[i]})" if nm.symmetry[i] else ""
            flag = "   [imaginary]" if nm.freq[i] < 0 else ""
            label = f"Mode {i + 1:3d}   {nm.freq[i]:9.1f} cm^-1{sym_label}{flag}"
            if order_name == "IR intensity":
                label += f"   IR={nm.IR[i]:.1f}"
            elif order_name == "Raman intensity":
                label += f"   Raman={nm.Raman[i]:.1f}"
            items.append(label)
        return items, [i + 1 for i in idx]

    # ---- build the selector window --------------------------------------
    # Reuse an already-existing Tk root/interpreter if one is active (e.g.
    # IPython/Spyder's own hidden Tk app that drives its "tk" event-loop
    # integration, or one matplotlib itself created) instead of calling
    # tk.Tk() again: a second independent Tcl/Tk interpreter in the same
    # process is a known source of native (segfault-level) crashes on
    # macOS, distinct from the earlier native-backend rotation crash.
    _existing_root = tk._default_root
    standalone_root = _existing_root is None
    root = tk.Tk() if standalone_root else tk.Toplevel(_existing_root)
    root.title(f"G16 Normal Mode Viewer - {filename}")

    screen_w = root.winfo_screenwidth()
    screen_h = root.winfo_screenheight()
    panel_w = 480
    root.geometry(f"{panel_w}x580+20+40")

    order_names = ["Mode number", "IR intensity"]
    if nm.has_Raman:
        order_names.append("Raman intensity")

    state["mode_items"], state["mode_data"] = build_mode_items("Mode number")

    frm = ttk.Frame(root, padding=10)
    frm.pack(fill="both", expand=True)

    ttk.Label(frm, text="Select a vibrational mode:", font=("", 10, "bold")).pack(anchor="w")
    mode_var = tk.StringVar(value=state["mode_items"][0])
    mode_combo = ttk.Combobox(frm, textvariable=mode_var, values=state["mode_items"],
                               state="readonly", width=55)
    mode_combo.pack(fill="x", pady=(0, 10))

    order_row = ttk.Frame(frm)
    order_row.pack(fill="x", pady=(0, 10))
    ttk.Label(order_row, text="Order by:").pack(side="left")
    order_var = tk.StringVar(value="Mode number")
    order_combo = ttk.Combobox(order_row, textvariable=order_var, values=order_names,
                                state="readonly", width=20)
    order_combo.pack(side="left", padx=(5, 0))

    ttk.Label(frm, text="Figure title:", font=("", 10, "bold")).pack(anchor="w")
    title_row = ttk.Frame(frm)
    title_row.pack(fill="x", pady=(0, 10))
    title_var = tk.StringVar()
    ttk.Entry(title_row, textvariable=title_var, width=40).pack(side="left", fill="x", expand=True)
    apply_title_btn = ttk.Button(title_row, text="Apply title")
    apply_title_btn.pack(side="left", padx=(5, 0))

    ttk.Label(frm, text="g16_draw_mode options:", font=("", 10, "bold")).pack(anchor="w", pady=(5, 0))

    opt1 = ttk.Frame(frm)
    opt1.pack(fill="x", pady=2)
    ttk.Label(opt1, text="Scale:", width=12).pack(side="left")
    scale_var = tk.StringVar(value=str(scale))
    scale_combo = ttk.Combobox(opt1, textvariable=scale_var, values=["0.5", "1", "1.5", "2", "3"], width=8)
    scale_combo.pack(side="left")
    ttk.Label(opt1, text="ArrowColor:", width=12).pack(side="left", padx=(15, 0))
    color_var = tk.StringVar(value="Custom..." if arrow_color is not None else "Orange (default)")
    color_combo = ttk.Combobox(opt1, textvariable=color_var, values=_COLOR_NAMES,
                                state="readonly", width=16)
    color_combo.pack(side="left")

    opt2 = ttk.Frame(frm)
    opt2.pack(fill="x", pady=2)
    ttk.Label(opt2, text="AtomScale:", width=12).pack(side="left")
    atom_scale_var = tk.StringVar(value=str(atom_scale))
    atom_scale_combo = ttk.Combobox(opt2, textvariable=atom_scale_var,
                                     values=["0.2", "0.3", "0.35", "0.5"], width=8)
    atom_scale_combo.pack(side="left")
    ttk.Label(opt2, text="BondTol:", width=12).pack(side="left", padx=(15, 0))
    bond_tol_var = tk.StringVar(value=str(bond_tol))
    bond_tol_combo = ttk.Combobox(opt2, textvariable=bond_tol_var,
                                   values=["1.10", "1.20", "1.30", "1.40", "1.50"], width=8)
    bond_tol_combo.pack(side="left")

    opt3 = ttk.Frame(frm)
    opt3.pack(fill="x", pady=2)
    show_labels_var = tk.BooleanVar(value=show_labels)
    show_labels_chk = ttk.Checkbutton(opt3, text="ShowLabels", variable=show_labels_var)
    show_labels_chk.pack(side="left")
    flip_sign_var = tk.BooleanVar(value=flip_sign)
    flip_sign_chk = ttk.Checkbutton(opt3, text="FlipSign", variable=flip_sign_var)
    flip_sign_chk.pack(side="left", padx=(15, 0))

    opt4 = ttk.Frame(frm)
    opt4.pack(fill="x", pady=2)
    ttk.Label(opt4, text="LabelFontSize:", width=12).pack(side="left")
    font_size_var = tk.StringVar(value=str(font_size))
    font_size_combo = ttk.Combobox(opt4, textvariable=font_size_var,
                                    values=["6", "8", "10", "12", "14", "16"], width=8)
    font_size_combo.pack(side="left")

    ttk.Label(frm, text="Save current mode figure as:", font=("", 10, "bold")).pack(anchor="w", pady=(10, 0))
    save_row = ttk.Frame(frm)
    save_row.pack(fill="x", pady=2)
    fmt_var = tk.StringVar(value="PDF (vector)")
    ttk.Combobox(save_row, textvariable=fmt_var, values=list(_FMT_TO_EXT),
                 state="readonly", width=18).pack(side="left")
    save_btn = ttk.Button(save_row, text="Save figure...")
    save_btn.pack(side="left", padx=(10, 0))

    ttk.Label(frm, text=f"{nm.Natoms} atoms  |  {nm.Nmodes} modes  |  file: {filename}",
              foreground="#737373", wraplength=panel_w - 20).pack(anchor="w", pady=(15, 0))

    # ---- callbacks ---------------------------------------------------------
    def get_selected_mode_idx():
        try:
            i = state["mode_items"].index(mode_var.get())
        except ValueError:
            i = 0
        return state["mode_data"][i]

    def build_draw_opts():
        color = state["custom_color"] if color_var.get() == "Custom..." else _COLOR_PRESETS[color_var.get()]
        opts = dict(
            scale=float(scale_var.get()),
            arrow_color=color,
            atom_scale=float(atom_scale_var.get()),
            bond_tol=float(bond_tol_var.get()),
            show_labels=bool(show_labels_var.get()),
            flip_sign=bool(flip_sign_var.get()),
        )
        opts.update(extra_opts)
        return opts

    def apply_label_font_size(ax, fsz):
        for t in ax.texts:
            t.set_fontsize(fsz)

    def position_mode_window(fig):
        cascade_step, cascade_max = 40, 6
        idx = len(state["open_figs"]) % cascade_max
        fig_w = max(500, screen_w - panel_w - 40 - cascade_max * cascade_step)
        fig_h = max(400, screen_h - 80 - cascade_max * cascade_step)
        fig_x = panel_w + 20 + idx * cascade_step
        fig_y = 40 + (cascade_max - idx) * cascade_step
        try:
            fig.canvas.manager.window.geometry(f"{fig_w}x{fig_h}+{fig_x}+{fig_y}")
        except Exception:
            pass

    def draw_mode(k, close_old):
        to_close = state["current_fig"] if (close_old and state["current_fig"] is not None) else None

        try:
            opts = build_draw_opts()
            ax = g16_draw_mode(mol, nm, k, **opts)
        except Exception as e:
            messagebox.showerror("g16_draw_mode error", str(e), parent=root)
            return

        fig = ax.figure
        try:
            apply_label_font_size(ax, float(font_size_var.get()))
        except ValueError:
            pass

        fig.canvas.manager.set_window_title(f"{filename} - Mode {k} ({nm.freq[k - 1]:.1f} cm-1)")
        position_mode_window(fig)
        fig.canvas.manager.show()

        if to_close is not None and to_close is not fig:
            import matplotlib.pyplot as plt
            plt.close(to_close)
            state["open_figs"] = [f for f in state["open_figs"] if f is not to_close]

        state["current_fig"] = fig
        state["open_figs"].append(fig)

    def on_mode_selected(event=None):
        draw_mode(get_selected_mode_idx(), close_old=False)

    def on_order_changed(event=None):
        cur_idx = get_selected_mode_idx()
        new_items, new_data = build_mode_items(order_var.get())
        state["mode_items"], state["mode_data"] = new_items, new_data
        mode_combo["values"] = new_items
        try:
            mode_var.set(new_items[new_data.index(cur_idx)])
        except ValueError:
            mode_var.set(new_items[0])

    def on_option_changed(event=None):
        draw_mode(get_selected_mode_idx(), close_old=True)

    def on_color_changed(event=None):
        if color_var.get() == "Custom...":
            rgb, _hex = colorchooser.askcolor(color=_to_hex(state["custom_color"]),
                                               title="Pick arrow colour", parent=root)
            if rgb is not None:
                state["custom_color"] = tuple(c / 255 for c in rgb)
        on_option_changed()

    def resolve_target_figure():
        return state["current_fig"]

    def apply_title(event=None):
        target = resolve_target_figure()
        if target is None:
            messagebox.showinfo("Nothing to title", "No mode figure is currently open.", parent=root)
            return
        new_title = title_var.get()
        if not new_title:
            return
        if target.axes:
            target.axes[0].set_title(new_title)
        target.canvas.manager.set_window_title(new_title)
        target.canvas.draw_idle()

    def save_current_figure():
        target = resolve_target_figure()
        if target is None:
            messagebox.showinfo("Nothing to save", "No mode figure is currently open.", parent=root)
            return
        ext = _FMT_TO_EXT[fmt_var.get()]
        out_file = filedialog.asksaveasfilename(
            title="Save mode figure as", defaultextension=ext,
            filetypes=[(fmt_var.get(), f"*{ext}")], initialfile=f"mode_figure{ext}", parent=root,
        )
        if not out_file:
            return
        try:
            target.savefig(out_file, dpi=300 if ext == ".jpg" else None)
            print(f"g16_mode_viewer: figure saved to {out_file}")
        except Exception as e:
            messagebox.showerror("Export error", str(e), parent=root)

    def on_close():
        import matplotlib.pyplot as plt
        for f in list(state["open_figs"]):
            plt.close(f)
        root.destroy()

    mode_combo.bind("<<ComboboxSelected>>", on_mode_selected)
    order_combo.bind("<<ComboboxSelected>>", on_order_changed)
    for w in (scale_combo, atom_scale_combo, bond_tol_combo, font_size_combo):
        w.bind("<<ComboboxSelected>>", on_option_changed)
        w.bind("<Return>", on_option_changed)
        w.bind("<FocusOut>", on_option_changed)
    color_combo.bind("<<ComboboxSelected>>", on_color_changed)
    show_labels_chk.config(command=on_option_changed)
    flip_sign_chk.config(command=on_option_changed)
    apply_title_btn.config(command=apply_title)
    save_btn.config(command=save_current_figure)
    root.protocol("WM_DELETE_WINDOW", on_close)

    draw_mode(state["mode_data"][0], close_old=False)

    if not _running_ipython_tk_loop():
        root.mainloop()
