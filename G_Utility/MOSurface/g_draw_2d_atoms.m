function g_draw_2d_atoms(ax, symbols, atom_xy, xyz_ang)
%G_DRAW_2D_ATOMS  Lightweight 2D sketch of atoms (CPK-coloured markers +
%   element/index labels) and simple (single-line, no bond-order)
%   distance-based bonds, projected onto a 'Mode','contour' cutting
%   plane. Shared by G_draw_mo_surface and G_draw_density_surface.
%
%   g_draw_2d_atoms(ax, symbols, atom_xy, xyz_ang)
%
%   ATOM_XY [Nat x 2] are the already-projected in-plane plotting
%   positions (Angstrom); XYZ_ANG [Nat x 3] are the true (unprojected) 3D
%   positions (Angstrom), used only for the bond-length distance test, so
%   a bond between two atoms that are close in 3D but far apart in the
%   projection is still drawn correctly.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

cpk_colors = containers.Map( ...
    {'H','C','N','O','F','P','S','Cl','Br','I','Au'}, ...
    {[0.60 0.80 1.00],[0.30 0.30 0.30],[0.10 0.30 0.90],[0.90 0.10 0.10], ...
     [0.20 0.80 0.20],[1.00 0.50 0.00],[1.00 0.85 0.00],[0.20 0.85 0.20], ...
     [0.55 0.20 0.10],[0.45 0.00 0.65],[1.00 0.82 0.14]});
cov_radii = containers.Map( ...
    {'H','C','N','O','F','P','S','Cl','Br','I','Au'}, ...
    {0.31,0.76,0.71,0.66,0.57,1.07,1.05,1.02,1.20,1.39,1.36});
default_color  = [0.65 0.20 0.80];
default_radius = 0.80;
bond_color     = [0.45 0.45 0.45];

Nat = numel(symbols);
for i = 1:Nat
    ri = local_radius(symbols{i});
    for j = i+1:Nat
        rj = local_radius(symbols{j});
        d = norm(xyz_ang(i,:) - xyz_ang(j,:));
        if d < (ri + rj) * 1.30
            plot(ax, [atom_xy(i,1) atom_xy(j,1)], [atom_xy(i,2) atom_xy(j,2)], ...
                '-', 'Color', bond_color, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
end

for i = 1:Nat
    clr = local_color(symbols{i});
    plot(ax, atom_xy(i,1), atom_xy(i,2), 'o', 'MarkerFaceColor', clr, ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'HandleVisibility', 'off');
    text(ax, atom_xy(i,1), atom_xy(i,2), sprintf('  %s%d', symbols{i}, i), ...
        'FontSize', 7, 'Color', clr*0.6, 'Interpreter', 'none', 'HandleVisibility', 'off');
end

    function clr = local_color(sym)
        if isKey(cpk_colors, sym), clr = cpk_colors(sym); else, clr = default_color; end
    end
    function r = local_radius(sym)
        if isKey(cov_radii, sym), r = cov_radii(sym); else, r = default_radius; end
    end

end % g_draw_2d_atoms
