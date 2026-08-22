import numpy as np

from .draw_molecule import g16_draw_molecule
from .get_bond_length import g16_get_bond_length


def g16_draw_deformation(mol1, mol2, overlay=False, scale=1, arrow_color=(1.0, 0.4, 0.1),
                          overlay2_color=(0.75, 0.1, 0.1), atom_scale=0.35, bond_tol=1.30,
                          show_labels=False, arrow_threshold=0.05):
    """Visualises the geometric deformation between two related structures
    of the same molecule (matplotlib quiver arrows in place of MATLAB's
    hand-built cone-tipped arrows -- same information, simpler static
    rendering, consistent with g16_draw_mode's own port).

    Compares two calculations of the same molecule -- typically with vs
    without a small static electric field (finite-field NLO workflows),
    two conformers, or before/after optimisation -- by drawing an arrow
    from each atom's position in mol1 to its position in mol2
    (mol2.xyz - mol1.xyz).

    Parameters
    ----------
    mol1, mol2 : Struct -- as returned by g16_structure (need .symbols,
        .xyz, .Natoms), same atom count and ordering.
    overlay : bool -- False (default): draw mol1 only, as a full CPK
        ball-and-stick model, with displacement arrows.
        True: ALSO draws mol2 as a simplified skeletal overlay (dashed
        lines + small dots, no spheres), in overlay2_color, drawn at
        mol1.xyz + scale*(mol2.xyz - mol1.xyz) -- the SAME exaggerated
        position the arrows point to, not mol2's true coordinates, so
        both ends of every arrow stay visually anchored. Bond
        connectivity for the overlay is still detected from mol2's true,
        unscaled geometry (an exaggerated geometry at a large scale is
        not a real molecular geometry and would misjudge bonding
        distances).
    scale : float -- arrow length AND overlay-position multiplier
        (default 1), applied identically to both, unlike g16_draw_mode's
        scale (which normalises to a fixed maximum arrow length): here
        the raw displacement is multiplied as-is, so the printed max
        displacement/RMSD tells you what scale is actually needed to see
        anything (real deformations are often only ~1e-3 to 1e-2
        Angstrom).
    arrow_color : tuple -- default (1.0, 0.4, 0.1), same as g16_draw_mode
    overlay2_color : tuple -- colour of mol2's skeletal overlay when
        overlay=True (default (0.75, 0.1, 0.1))
    atom_scale, bond_tol, show_labels : see g16_draw_molecule (applied to
        mol1's full CPK render)
    arrow_threshold : float -- skip the arrow for atoms whose
        displacement is below this fraction of the single largest
        per-atom displacement (default 0.05)

    Prints to the console: the largest single-atom displacement (with
    which atom), and the RMSD over all atoms -- the numbers to look at
    before picking `scale`.

    Returns
    -------
    ax : the matplotlib 3D axes used

    Example
    -------
    mol1 = g16_structure('nofield.out')
    mol2 = g16_structure('field_x025.out')
    g16_draw_deformation(mol1, mol2, scale=200)
    g16_draw_deformation(mol1, mol2, overlay=True, scale=200)

    See also g16_draw_molecule, g16_draw_mode.
    """
    if not hasattr(mol1, "xyz") or not hasattr(mol2, "xyz"):
        raise ValueError("g16_draw_deformation: mol1 and mol2 must both be mol structs with .xyz.")
    if mol1.Natoms != mol2.Natoms:
        raise ValueError(
            f"g16_draw_deformation: mol1.Natoms ({mol1.Natoms}) does not match "
            f"mol2.Natoms ({mol2.Natoms}) -- not the same molecule."
        )

    xyz1 = np.asarray(mol1.xyz, dtype=float)
    xyz2 = np.asarray(mol2.xyz, dtype=float)

    U = xyz2 - xyz1  # [Natoms x 3], Angstrom
    norms_i = np.linalg.norm(U, axis=1)
    i_max = int(np.argmax(norms_i))
    max_d = float(norms_i[i_max])
    rmsd = float(np.sqrt(np.mean(np.sum(U**2, axis=1))))

    print("\n-- g16_draw_deformation --")
    if getattr(mol1, "symbols", None) is not None and len(mol1.symbols) > i_max:
        print(f"  Max displacement : {max_d:.5f} A  (atom {i_max + 1}, {mol1.symbols[i_max]})")
    else:
        print(f"  Max displacement : {max_d:.5f} A  (atom {i_max + 1})")
    print(f"  RMSD             : {rmsd:.5f} A\n")

    if max_d == 0:
        import warnings
        warnings.warn("g16_draw_deformation: mol1 and mol2 have identical coordinates -- nothing to draw.")
        return None

    import matplotlib.pyplot as plt

    fig = plt.figure()
    fig.canvas.manager.set_window_title("Structural deformation")
    ax = fig.add_subplot(111, projection="3d")

    title_str = f"Deformation  --  max {max_d:.4f} A, RMSD {rmsd:.4f} A"
    g16_draw_molecule(mol1, ax=ax, atom_scale=atom_scale, bond_tol=bond_tol,
                       show_labels=show_labels, show_legend=False, title=title_str)

    # -------------------------------------------------------------------
    # Optional skeletal overlay of mol2, drawn at the SAME exaggerated
    # position as the arrows below (mol1.xyz + scale*U) -- see docstring.
    # -------------------------------------------------------------------
    U_scaled = U * scale
    xyz2_draw = xyz1 + U_scaled

    if overlay:
        bond_table2 = g16_get_bond_length(mol2, tolerance=bond_tol, include_h=True)
        for _, row in bond_table2.iterrows():
            a1 = int(row["Atom1"]) - 1
            a2 = int(row["Atom2"]) - 1
            ax.plot(xyz2_draw[[a1, a2], 0], xyz2_draw[[a1, a2], 1], xyz2_draw[[a1, a2], 2],
                    color=overlay2_color, linewidth=1.2, linestyle="--")
        ax.scatter(xyz2_draw[:, 0], xyz2_draw[:, 1], xyz2_draw[:, 2],
                   s=16, color=overlay2_color, depthshade=False)

    # -------------------------------------------------------------------
    # Displacement arrows, mol1(i) -> mol2(i)
    # -------------------------------------------------------------------
    for i in range(mol1.Natoms):
        if norms_i[i] / max_d < arrow_threshold:
            continue
        x0, y0, z0 = xyz1[i]
        dx, dy, dz = U_scaled[i]
        ax.quiver(x0, y0, z0, dx, dy, dz, color=arrow_color, linewidth=2.0,
                  arrow_length_ratio=0.25)

    return ax
