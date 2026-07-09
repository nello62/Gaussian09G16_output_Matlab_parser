import os

from ._common import HARTREE_TO_EV


def g16_draw_orbital(oe, nlevels=5, units="eV", ax=None,
                      occ_color=(0.25, 0.25, 0.70), virt_color=(0.70, 0.30, 0.25),
                      homo_color=(0.00, 0.45, 0.85), lumo_color=(0.90, 0.35, 0.00),
                      arrow_color=(0.15, 0.60, 0.15), show_labels=True, show_spins=True,
                      font_size=8, frontier_font_size=13, arrow_head_size=0.2, title=""):
    """Draws a molecular-orbital energy-level diagram from the Struct
    returned by g16_orbital_energies, highlighting the HOMO-LUMO
    transition with an arrow labelled with the gap.

    Parameters mirror G16_draw_orbital.m — see its docstring for details.

    Returns
    -------
    ax : the matplotlib axes used
    """
    import matplotlib.pyplot as plt

    if not hasattr(oe, "alpha_occ") or not hasattr(oe, "alpha_virt"):
        raise ValueError("g16_draw_orbital: oe must be the Struct returned by g16_orbital_energies.")
    if oe.alpha_occ.size == 0 or oe.alpha_virt.size == 0:
        raise ValueError("g16_draw_orbital: oe has no occupied/virtual orbital energies to plot.")

    units_key = units.lower()
    if units_key == "ev":
        occ_e = oe.alpha_occ * HARTREE_TO_EV
        virt_e = oe.alpha_virt * HARTREE_TO_EV
        gap_disp = oe.gap_eV
        unit_lbl = "eV"
    elif units_key == "hartree":
        occ_e = oe.alpha_occ
        virt_e = oe.alpha_virt
        gap_disp = oe.gap
        unit_lbl = "Ha"
    else:
        raise ValueError("g16_draw_orbital: units must be 'eV' or 'Hartree'.")

    nocc, nvirt = len(occ_e), len(virt_e)
    occ_show = occ_e[max(0, nocc - nlevels):nocc]
    virt_show = virt_e[:min(nlevels, nvirt)]

    HOMO_e = occ_e[-1]
    LUMO_e = virt_e[0]

    if ax is None:
        fig, ax = plt.subplots()
        fig.canvas.manager.set_window_title("Orbital energy diagram")

    bar_x = (-0.5, 0.5)

    for i, e in enumerate(occ_show):
        is_homo = i == len(occ_show) - 1
        clr, lw = (homo_color, 2.5) if is_homo else (occ_color, 1.5)
        ax.plot(bar_x, [e, e], color=clr, linewidth=lw)
        if show_spins:
            ax.text(0, e, r"$\uparrow\downarrow$", fontsize=font_size + 2, color=clr,
                    ha="center", va="center")
        if show_labels:
            ax.text(bar_x[0] - 0.08, e, f"{e:.3f}", fontsize=font_size, color=clr,
                    ha="right", va="center")

    for i, e in enumerate(virt_show):
        is_lumo = i == 0
        clr, lw = (lumo_color, 2.5) if is_lumo else (virt_color, 1.5)
        ax.plot(bar_x, [e, e], color=clr, linewidth=lw)
        if show_labels:
            ax.text(bar_x[0] - 0.08, e, f"{e:.3f}", fontsize=font_size, color=clr,
                    ha="right", va="center")

    ax.text(bar_x[1] + 0.08, HOMO_e, "HOMO", color=homo_color, fontweight="bold",
            fontsize=frontier_font_size, ha="left", va="center")
    ax.text(bar_x[1] + 0.08, LUMO_e, "LUMO", color=lumo_color, fontweight="bold",
            fontsize=frontier_font_size, ha="left", va="center")

    arrow_x = 0.95
    ax.annotate("", xy=(arrow_x, LUMO_e), xytext=(arrow_x, HOMO_e),
                arrowprops=dict(color=arrow_color, linewidth=2,
                                headwidth=8 * arrow_head_size / 0.2,
                                headlength=10 * arrow_head_size / 0.2))
    ax.text(arrow_x + 0.12, (HOMO_e + LUMO_e) / 2, f"$\\Delta E$ = {gap_disp:.3f} {unit_lbl}",
            fontsize=frontier_font_size, fontweight="bold", color=arrow_color,
            ha="left", va="center")

    ax.set_xlim(-1.0, 2.3)
    ax.set_xticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.set_ylabel(f"Energy ({unit_lbl})")

    if not title:
        if getattr(oe, "filename", None):
            title = os.path.splitext(os.path.basename(oe.filename))[0]
        else:
            title = "Orbital Energy Diagram"
    ax.set_title(title, fontsize=11)

    return ax
