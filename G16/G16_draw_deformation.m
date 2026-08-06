function G16_draw_deformation(mol1, mol2, varargin)
% G16_DRAW_DEFORMATION  Visualises the geometric deformation between two
%                       related structures of the same molecule.
%
%   G16_draw_deformation(mol1, mol2)
%   G16_draw_deformation(mol1, mol2, 'Name', Value, ...)
%
%   mol1/mol2 - structs with .symbols/.xyz/.Natoms (as returned by
%               G16_structure, G_read_structure_file, G16_fchk_read's
%               .mol sub-struct, ...), same atom count and ordering.
%
%   Compares two calculations of the same molecule -- typically with vs
%   without a small static electric field (finite-field NLO workflows),
%   two conformers, or before/after optimisation -- by drawing an arrow
%   from each atom's position in mol1 to its position in mol2
%   (mol2.xyz - mol1.xyz), the same displacement-arrow style used by
%   G16_DRAW_MODE for a normal mode, but here the actual geometric shift
%   between two structures rather than a vibrational eigenvector.
%
%   Optional parameters:
%       'Overlay'        - false (default): draw mol1 only, as a full CPK
%                          ball-and-stick model, with displacement arrows.
%                        - true: ALSO draws mol2 as a simplified skeletal
%                          overlay (thin lines + small dots, no spheres),
%                          in 'Overlay2Color', so both ends of every arrow
%                          are visually anchored to a real atom position
%                          instead of just an arrowhead in empty space.
%                          Drawn at mol1.xyz + Scale*(mol2.xyz-mol1.xyz),
%                          i.e. the SAME exaggerated position the arrows
%                          point to, not mol2's true coordinates -- so the
%                          overlay always lines up with the arrow tips.
%                          Bond connectivity for the overlay is still
%                          detected from mol2's true, unscaled geometry.
%       'Scale'          - arrow length AND overlay-position multiplier
%                          (default 1), applied identically to both so
%                          they always stay visually consistent. Real
%                          deformations, e.g. field-induced geometry
%                          shifts, are often only ~1e-3 to 1e-2 Angstrom --
%                          far too small to see at Scale=1. The console
%                          report (see below) prints the true maximum
%                          displacement so you can pick a sensible value.
%       'ArrowColor'     - (default [1.0 0.4 0.1], same as G16_draw_mode)
%       'Overlay2Color'  - colour of mol2's skeletal overlay when
%                          'Overlay' is true (default [0.75 0.1 0.1])
%       'AtomScale'      - CPK sphere scale for mol1 (default 0.35)
%       'BondTol'        - bond-detection tolerance, both structures
%                          (default 1.30)
%       'ShowLabels'     - show atom index labels on mol1 (default false)
%       'ArrowThreshold' - skip the arrow for atoms whose displacement is
%                          below this fraction of the single largest
%                          per-atom displacement (default 0.05)
%
%   Prints to the console: the largest single-atom displacement (with
%   which atom), and the RMSD over all atoms -- the numbers to look at
%   before picking 'Scale'.
%
%   Example:
%       mol1 = G16_structure('nofield.out');
%       mol2 = G16_structure('field_x025.out');
%       G16_draw_deformation(mol1, mol2, 'Scale', 200);
%       G16_draw_deformation(mol1, mol2, 'Overlay', true, 'Scale', 200);
%
%   See also G16_DRAW_MOLECULE, G16_DRAW_MODE, G_MATCH_MODES.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

p = inputParser;
addRequired(p,  'mol1');
addRequired(p,  'mol2');
addParameter(p, 'Overlay',        false,          @islogical);
addParameter(p, 'Scale',          1,              @isnumeric);
addParameter(p, 'ArrowColor',     [1.0 0.4 0.1],  @isnumeric);
addParameter(p, 'Overlay2Color',  [0.75 0.1 0.1], @isnumeric);
addParameter(p, 'AtomScale',      0.35,           @isnumeric);
addParameter(p, 'BondTol',        1.30,           @isnumeric);
addParameter(p, 'ShowLabels',     false,          @islogical);
addParameter(p, 'ArrowThreshold', 0.05,           @isnumeric);
parse(p, mol1, mol2, varargin{:});

do_overlay   = p.Results.Overlay;
scale        = p.Results.Scale;
arrow_color  = p.Results.ArrowColor;
ov2_color    = p.Results.Overlay2Color;
atom_scale   = p.Results.AtomScale;
bond_tol     = p.Results.BondTol;
show_labels  = p.Results.ShowLabels;
thresh_frac  = p.Results.ArrowThreshold;

if ~isstruct(mol1) || ~isfield(mol1, 'xyz') || ~isstruct(mol2) || ~isfield(mol2, 'xyz')
    error('G16_draw_deformation: mol1 and mol2 must both be mol structs with .xyz.');
end
if mol1.Natoms ~= mol2.Natoms
    error('G16_draw_deformation: mol1.Natoms (%d) does not match mol2.Natoms (%d) -- not the same molecule.', ...
        mol1.Natoms, mol2.Natoms);
end

% -------------------------------------------------------------------------
% Displacement and console summary
% -------------------------------------------------------------------------
U = mol2.xyz - mol1.xyz;                  % [Natoms x 3], Angstrom
norms_i = sqrt(sum(U.^2, 2));
[max_d, i_max] = max(norms_i);
rmsd = sqrt(mean(sum(U.^2, 2)));

fprintf('\n── G16_draw_deformation ──\n');
if isfield(mol1, 'symbols') && numel(mol1.symbols) >= i_max
    fprintf('  Max displacement : %.5f A  (atom %d, %s)\n', max_d, i_max, mol1.symbols{i_max});
else
    fprintf('  Max displacement : %.5f A  (atom %d)\n', max_d, i_max);
end
fprintf('  RMSD             : %.5f A\n\n', rmsd);

if max_d == 0
    warning('G16_draw_deformation: mol1 and mol2 have identical coordinates -- nothing to draw.');
    return
end

% -------------------------------------------------------------------------
% Base render: mol1 as a full CPK structure
% -------------------------------------------------------------------------
fig = figure('Color', 'white', 'Name', 'Structural deformation', 'NumberTitle', 'off');
ax = axes('Parent', fig);

title_str = sprintf('Deformation  --  max %.4f A, RMSD %.4f A', max_d, rmsd);
G16_draw_molecule(mol1, 'Ax', ax, 'AtomScale', atom_scale, 'BondTol', bond_tol, ...
    'ShowLabels', show_labels, 'ShowLegend', false, 'Title', title_str);
hold(ax, 'on');

% -------------------------------------------------------------------------
% Optional skeletal overlay of mol2 (thin lines + small dots, no spheres)
% -------------------------------------------------------------------------
% Drawn at the SAME exaggerated position as the arrows below
% (mol1.xyz + scale*U), not at mol2's true, unscaled coordinates: with
% 'Scale' ~= 1 (routine for these typically-tiny deformations), drawing
% the overlay at the true position would leave it visually detached from
% the arrow tips, which point at the exaggerated position instead.
% Connectivity (which atoms count as bonded) is still detected from the
% true, unscaled mol2.xyz, since an exaggerated geometry at a large Scale
% is not a real molecular geometry and would misjudge bonding distances.
% -------------------------------------------------------------------------
U_scaled = U * scale;
xyz2_draw = mol1.xyz + U_scaled;

if do_overlay
    bondTable2 = G16_get_bond_length(mol2, 'Tolerance', bond_tol, 'IncludeH', true);
    for b = 1:height(bondTable2)
        a1 = bondTable2.Atom1(b);
        a2 = bondTable2.Atom2(b);
        line(ax, xyz2_draw([a1 a2], 1), xyz2_draw([a1 a2], 2), xyz2_draw([a1 a2], 3), ...
            'Color', ov2_color, 'LineWidth', 1.2, 'LineStyle', '--', 'HandleVisibility', 'off');
    end
    plot3(ax, xyz2_draw(:,1), xyz2_draw(:,2), xyz2_draw(:,3), ...
        'o', 'MarkerSize', 4, 'MarkerFaceColor', ov2_color, 'MarkerEdgeColor', ov2_color, ...
        'HandleVisibility', 'off');
end

% -------------------------------------------------------------------------
% Displacement arrows, mol1(i) -> mol2(i)
% -------------------------------------------------------------------------
for i = 1:mol1.Natoms
    if norms_i(i) / max_d < thresh_frac
        continue
    end
    draw_arrow3(ax, mol1.xyz(i,1), mol1.xyz(i,2), mol1.xyz(i,3), ...
        U_scaled(i,1), U_scaled(i,2), U_scaled(i,3), arrow_color);
end

rotate3d(ax, 'on');

end % G16_draw_deformation


% =========================================================================
%  3D arrow with cone tip (same construction as G16_draw_mode's own
%  local helper, duplicated here since it is not a shared/exported
%  function)
% =========================================================================
function draw_arrow3(ax, x0, y0, z0, dx, dy, dz, color)
line(ax, [x0, x0+dx], [y0, y0+dy], [z0, z0+dz], ...
     'Color', color, 'LineWidth', 2.0, 'HandleVisibility', 'off');

len = sqrt(dx^2 + dy^2 + dz^2);
if len < 1e-6, return; end

tip_frac = 0.25;
tip_r    = 0.07;

tip_len = len * tip_frac;
ux = dx/len; uy = dy/len; uz = dz/len;

cx = x0 + dx - ux*tip_len;
cy = y0 + dy - uy*tip_len;
cz = z0 + dz - uz*tip_len;

if abs(ux) < 0.9
    perp = [0 -uz uy];
else
    perp = [-uz 0 ux];
end
perp = perp / norm(perp);
perp2 = cross([ux uy uz], perp);

th = linspace(0, 2*pi, 16);
bx = cx + tip_r*(cos(th)*perp(1) + sin(th)*perp2(1));
by = cy + tip_r*(cos(th)*perp(2) + sin(th)*perp2(2));
bz = cz + tip_r*(cos(th)*perp(3) + sin(th)*perp2(3));

tip_x = x0 + dx;
tip_y = y0 + dy;
tip_z = z0 + dz;

for j = 1:numel(th)-1
    patch(ax, ...
        [bx(j), bx(j+1), tip_x], ...
        [by(j), by(j+1), tip_y], ...
        [bz(j), bz(j+1), tip_z], ...
        color, 'EdgeColor', 'none', ...
        'FaceLighting', 'gouraud', ...
        'HandleVisibility', 'off');
end
end
