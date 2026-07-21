function G09_draw_molecule(mol, varargin)
% G09_DRAW_MOLECULE  Renders a 3D CPK ball-and-stick model from the mol structstructure.
%
%   G09_draw_molecule(mol)
%   G09_draw_molecule(mol, Name, Value, ...)
%
%   Input:
%       mol  - struct returned by G09_structure, with fields:
%                .symbols  {Natoms×1 cell}  simboli atomici
%                .xyz      [Natoms×3]       coordinate in Angstrom
%                .Natoms   int
%
%   Optional parameters (Name-Value):
%       'AtomScale'    - sphere radius as fraction of covalent radius (default: 0.35)
%       'BondTol'      - bond tolerance: bonded if dist < (r1+r2)*BondTol (default: 1.30)
%       'ShowLabels'   - show atom labels "C1", "H2", ... (default: true)
%       'ShowLegend'   - show element legend (default: true)
%       'Title'        - figure title (default: filename or 'Molecule')
%       'BgColor'      - axes background colour (default: [0.95 0.95 0.95])
%       'Ax'           - existing axes handle (default: nuova figura)
%       'ShowAxes'     - draw a small Cartesian X/Y/Z reference triad,
%                        anchored at the lower-left-front corner of the
%                        plotted molecule (default: false)
%       'AxesLength'   - length of the X/Y/Z arrows, in Angstrom (default:
%                        [] = auto, 20% of the plotted molecule's
%                        bounding-box diagonal)
%       'BondList'     - [Nbonds x 2] or [Nbonds x 3] explicit atom-index
%                        pairs (optionally with a 3rd column giving a
%                        pre-computed bond order 1/2/3) to draw as bonds,
%                        bypassing the BondTol distance criterion
%                        (default: [] = auto-detect from BondTol). Useful
%                        to keep a fixed bond topology (and, with the 3rd
%                        column, a fixed bond order) across a series of
%                        frames where atoms move (e.g. G09_ANIMATE_MODE),
%                        so bonds do not appear/disappear or flicker
%                        between single/double/triple as instantaneous
%                        distances change.
%
%   Bond order (single/double/triple) is estimated purely from bond
%   length for C-C, C-N, and C-O pairs (any other element pair is always
%   drawn as a single bond), and rendered as 1/2/3 parallel lines in the
%   usual chemical-drawing convention. This is a geometric estimate, not
%   Gaussian's own bond-order analysis (e.g. Wiberg/NBO indices).
%
%   Example:
%       mol = G09_structure('zeatin.out');
%       G09_draw_molecule(mol)
%       G09_draw_molecule(mol, 'ShowLabels', false, 'AtomScale', 0.4)
%       G09_draw_molecule(mol, 'ShowAxes', true)

% -------------------------------------------------------------------------
% Parse input arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'mol');
addParameter(p, 'AtomScale',  0.35,              @isnumeric);
addParameter(p, 'BondTol',    1.30,              @isnumeric);
addParameter(p, 'ShowLabels', true,              @islogical);
addParameter(p, 'ShowLegend', true,              @islogical);
addParameter(p, 'Title',      '',                @ischar);
addParameter(p, 'BgColor',    [0.95 0.95 0.95], @isnumeric);
addParameter(p, 'Ax',         [],                @ishandle);
addParameter(p, 'ShowAxes',   false,             @islogical);
addParameter(p, 'AxesLength', [],                @isnumeric);
addParameter(p, 'BondList',   [],                @isnumeric);
parse(p, mol, varargin{:});

atom_scale   = p.Results.AtomScale;
bond_tol     = p.Results.BondTol;
show_labels  = p.Results.ShowLabels;
show_legend  = p.Results.ShowLegend;
fig_title    = p.Results.Title;
bg_color     = p.Results.BgColor;
ax           = p.Results.Ax;
show_axes    = p.Results.ShowAxes;
axes_length  = p.Results.AxesLength;
bond_list    = p.Results.BondList;

% Validate mol struct
if ~isstruct(mol) || ~isfield(mol,'symbols') || ~isfield(mol,'xyz')
    error('G09_draw_molecule: mol must be the struct returned by G09_structure.');
end

% Default title
if isempty(fig_title)
    if isfield(mol, 'filename')
        [~, fn] = fileparts(mol.filename);
        fig_title = strrep(fn, '_', '\_');
    else
        fig_title = 'Molecule';
    end
end

% -------------------------------------------------------------------------
% CPK colour and radius tables
% -------------------------------------------------------------------------
cpk_colors = containers.Map( ...
    {'H',  'C',  'N',  'O',  'F',  'P',  'S',  'Cl', 'Br', 'I', ...
     'Au', 'Ag', 'Fe', 'Zn', 'Ca', 'Mg', 'Na', 'K',  'Si', 'B', 'Cu'}, ...
    {[0.60 0.80 1.00], ...   % H  - azzurro chiaro
     [0.30 0.30 0.30], ...   % C  - grigio scuro
     [0.10 0.30 0.90], ...   % N  - blu
     [0.90 0.10 0.10], ...   % O  - rosso
     [0.20 0.80 0.20], ...   % F  - verde
     [1.00 0.50 0.00], ...   % P  - arancio
     [1.00 0.85 0.00], ...   % S  - giallo
     [0.20 0.85 0.20], ...   % Cl - verde chiaro
     [0.55 0.20 0.10], ...   % Br - marrone rossiccio
     [0.45 0.00 0.65], ...   % I  - viola
     [1.00 0.82 0.14], ...   % Au - oro
     [0.75 0.75 0.75], ...   % Ag - argento
     [0.80 0.40 0.00], ...   % Fe - arancio-marrone
     [0.50 0.70 0.50], ...   % Zn - verde acqua
     [0.60 0.60 0.60], ...   % Ca - grigio
     [0.50 0.80 0.20], ...   % Mg - verde lime
     [0.65 0.40 0.90], ...   % Na - viola chiaro
     [0.55 0.20 0.85], ...   % K  - viola
     [0.60 0.50 0.40], ...   % Si - beige
     [1.00 0.65 0.50], ...   % B  - salmone
     [0.72 0.45 0.20]});     % Cu - rame

cov_radii = containers.Map( ...
    {'H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I', ...
     'Au','Ag','Fe','Zn','Ca','Mg','Na','K', 'Si','B', 'Cu'}, ...
    {0.31, 0.76, 0.71, 0.66, 0.57, 1.07, 1.05, 1.02, 1.20, 1.39, ...
     1.36, 1.45, 1.32, 1.22, 1.76, 1.41, 1.66, 2.03, 1.11, 0.84, 1.32});

default_color  = [0.65 0.20 0.80];  % purple for elements not in the CPK table
default_radius = 0.80;

% -------------------------------------------------------------------------
% Set up figure and axes
% -------------------------------------------------------------------------
if isempty(ax)
    fig = figure('Color', 'white', 'Name', fig_title, ...
                 'NumberTitle', 'off');
    ax = axes('Parent', fig);
end

hold(ax, 'on');
axis(ax, 'equal');
axis(ax, 'off');
set(ax, 'Color', bg_color);
view(ax, 3);
lighting(ax, 'gouraud');
material(ax, 'dull');
camlight(ax, 'headlight');
camlight(ax, 45, 30);

% -------------------------------------------------------------------------
% Draw bonds
% -------------------------------------------------------------------------
[th, ph] = meshgrid(linspace(0,2*pi,12), linspace(0,pi,8));  % low-resolution spheres for bonds

bond_color = [0.50 0.50 0.50];
if isempty(bond_list)
    for i = 1 : mol.Natoms
        ri = get_radius_local(mol.symbols{i}, cov_radii, default_radius);
        for j = i+1 : mol.Natoms
            rj = get_radius_local(mol.symbols{j}, cov_radii, default_radius);
            d  = norm(mol.xyz(i,:) - mol.xyz(j,:));
            if d < (ri + rj) * bond_tol
                order = classify_bond_order(mol.symbols{i}, mol.symbols{j}, d);
                draw_bond_lines(ax, mol.xyz(i,:), mol.xyz(j,:), order, bond_color);
            end
        end
    end
else
    for b = 1 : size(bond_list, 1)
        i = bond_list(b, 1);
        j = bond_list(b, 2);
        if size(bond_list, 2) >= 3
            order = bond_list(b, 3);   % pre-computed (e.g. by G09_animate_mode, to
                                        % keep the order fixed across frames)
        else
            d = norm(mol.xyz(i,:) - mol.xyz(j,:));
            order = classify_bond_order(mol.symbols{i}, mol.symbols{j}, d);
        end
        draw_bond_lines(ax, mol.xyz(i,:), mol.xyz(j,:), order, bond_color);
    end
end

% -------------------------------------------------------------------------
% Draw atoms (spheres) and collect legend handles
% -------------------------------------------------------------------------
[xs, ys, zs] = sphere(30);   % sfera ad alta risoluzione

syms_present = unique(mol.symbols, 'stable');
% sort: heavy elements first, H last
heavy = syms_present(~strcmp(syms_present,'H'));
if any(strcmp(syms_present,'H'))
    syms_ordered = [heavy; {'H'}];
else
    syms_ordered = heavy;
end

legend_handles = gobjects(numel(syms_ordered), 1);

for ki = 1 : numel(syms_ordered)
    sym = syms_ordered{ki};
    clr = get_color_local(sym, cpk_colors, default_color);
    r   = get_radius_local(sym, cov_radii, default_radius) * atom_scale;

    idx_atoms = find(strcmp(mol.symbols, sym));
    h_first = [];

    for ii = 1 : numel(idx_atoms)
        i   = idx_atoms(ii);
        cx  = mol.xyz(i,1);
        cy  = mol.xyz(i,2);
        cz  = mol.xyz(i,3);

        if ii == 1
            % First atom of this element: HandleVisibility on for legend
            h = surf(ax, xs*r + cx, ys*r + cy, zs*r + cz, ...
                     'FaceColor', clr, 'EdgeColor', 'none', ...
                     'FaceLighting', 'gouraud', ...
                     'DisplayName', sym);
            legend_handles(ki) = h;
        else
            % Subsequent atoms: do not appear in legend
            surf(ax, xs*r + cx, ys*r + cy, zs*r + cz, ...
                 'FaceColor', clr, 'EdgeColor', 'none', ...
                 'FaceLighting', 'gouraud', ...
                 'HandleVisibility', 'off');
        end

        % Etichetta
        if show_labels
            lbl_offset = r * 1.4;
            text(ax, cx + lbl_offset, cy + lbl_offset, cz + lbl_offset, ...
                 sprintf('%s%d', sym, i), ...
                 'FontSize', 7, 'Color', clr * 0.7, ...
                 'HorizontalAlignment', 'left', ...
                 'Interpreter', 'none', ...
                 'HandleVisibility', 'off');
        end
    end
end

% -------------------------------------------------------------------------
% Legend and title
% -------------------------------------------------------------------------
if show_legend
    leg = legend(ax, legend_handles, syms_ordered, ...
                 'Location', 'northwest', ...
                 'FontSize',  9, ...
                 'FontWeight','bold', ...
                 'Box',       'off', ...
                 'TextColor', [0.15 0.15 0.15]);
end

title(ax, fig_title, 'Interpreter', 'tex', 'FontSize', 11);

% Fit axis limits to the actual plotted molecule before placing the
% Cartesian axes indicator, so its corner anchor and default length are
% based on the true displayed extent.
axis(ax, 'tight');

% -------------------------------------------------------------------------
% Cartesian reference axes (optional) — small X/Y/Z indicator anchored at
% the lower-left-front corner of the plotted molecule (min corner of the
% axes limits), out of the way of the structure itself.
% -------------------------------------------------------------------------
if show_axes
    corner = [ax.XLim(1), ax.YLim(1), ax.ZLim(1)];
    diag_len = norm([diff(ax.XLim), diff(ax.YLim), diff(ax.ZLim)]);
    if isempty(axes_length)
        axes_length = diag_len * 0.20;
        if axes_length == 0 || isnan(axes_length)
            axes_length = 2;   % fallback, e.g. a single atom at the origin
        end
    end
    draw_cartesian_axes(ax, corner, axes_length);
end

% Enable interactive rotation
rotate3d(ax, 'on');

% Adjust perspective
camproj(ax, 'perspective');

end % function G09_draw_molecule


% =========================================================================
%  Local functions
% =========================================================================

function clr = get_color_local(sym, cpk_colors, default_color)
    if isKey(cpk_colors, sym)
        clr = cpk_colors(sym);
    else
        clr = default_color;
    end
end

function r = get_radius_local(sym, cov_radii, default_radius)
    if isKey(cov_radii, sym)
        r = cov_radii(sym);
    else
        r = default_radius;
    end
end

function order = classify_bond_order(sym_i, sym_j, d)
%CLASSIFY_BOND_ORDER  Estimates bond order (1/2/3) from bond length alone,
%   for C-C, C-N and C-O pairs; any other element pair is always treated
%   as a single bond. Purely geometric, like the rest of this toolbox's
%   bond-detection logic -- not derived from an actual Gaussian bond-order
%   analysis (e.g. Wiberg/NBO indices).
%
%   Thresholds are the midpoint between adjacent reference bond lengths
%   (triple/double/single, in Angstrom), EXCEPT the C-C double/single
%   boundary, which is set to 1.36 -- deliberately below the ~1.39-1.40 A
%   aromatic C-C range (verified on real ring systems), so symmetric
%   aromatic rings are drawn as all-single rather than all-double: real
%   aromatic bonds have no length alternation to recover from geometry
%   alone (bond order really is ~1.5 all around the ring), so "all
%   single" is the more honest rendering than "all double".
    pair = sort_pair_local(upper(sym_i), upper(sym_j));
    switch pair
        case 'CC'
            thresh = [1.27, 1.36];    % [triple/double, double/single]
        case 'CN'
            thresh = [1.22, 1.375];
        case 'CO'
            thresh = [1.165, 1.315];
        otherwise
            order = 1;
            return
    end
    if d < thresh(1)
        order = 3;
    elseif d < thresh(2)
        order = 2;
    else
        order = 1;
    end
end

function key = sort_pair_local(a, b)
%SORT_PAIR_LOCAL  Canonical (order-independent) 2-letter key for an
%   element pair, e.g. ('N','C') and ('C','N') both give 'CN'.
    p = sort({a, b});
    key = [p{1}, p{2}];
end

function draw_bond_lines(ax, p1, p2, order, color)
%DRAW_BOND_LINES  Draws ORDER (1, 2, or 3) parallel line segments between
%   P1 and P2, offset perpendicular to the bond axis, in the classic
%   double-/triple-bond drawing convention.
    u = p2 - p1;
    ulen = norm(u);
    if ulen == 0
        return
    end
    u = u / ulen;
    ref = [0 0 1];
    if abs(dot(u, ref)) > 0.9
        ref = [0 1 0];
    end
    v = cross(u, ref);
    v = v / norm(v);

    offset = 0.09;   % Angstrom, spacing between parallel bond lines
    switch order
        case 2
            shifts = [-1, 1] * (offset / 2);
        case 3
            shifts = [-1, 0, 1] * offset;
        otherwise
            shifts = 0;
    end

    for k = 1:numel(shifts)
        dv = v * shifts(k);
        line(ax, [p1(1) p2(1)] + dv(1), [p1(2) p2(2)] + dv(2), [p1(3) p2(3)] + dv(3), ...
             'Color', color, 'LineWidth', 2.0, 'HandleVisibility', 'off');
    end
end

function draw_cartesian_axes(ax, origin, axes_length)
%DRAW_CARTESIAN_AXES  Draws X/Y/Z reference arrows anchored at ORIGIN
%   (classic red/green/blue convention), each labelled at its tip.
    dirs   = eye(3) * axes_length;
    colors = {[0.85 0.10 0.10], [0.10 0.65 0.10], [0.10 0.10 0.85]};
    labels = {'X', 'Y', 'Z'};
    for a = 1:3
        d = dirs(a, :);
        quiver3(ax, origin(1), origin(2), origin(3), d(1), d(2), d(3), 0, ...
            'Color', colors{a}, 'LineWidth', 1.5, 'MaxHeadSize', 0.5, ...
            'HandleVisibility', 'off');
        text(ax, origin(1)+d(1), origin(2)+d(2), origin(3)+d(3), ...
             sprintf('  %s', labels{a}), ...
             'Color', colors{a}, 'FontSize', 10, 'FontWeight', 'bold', ...
             'Interpreter', 'none', 'HandleVisibility', 'off');
    end
end
