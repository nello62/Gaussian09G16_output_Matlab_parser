function G16_draw_mode(mol, nm, mode_idx, varargin)
% G16_DRAW_MODE  Visualises a normal mode on a 3D molecular structure.
%
%   G16_draw_mode(mol, nm, mode_idx)
%   G16_draw_mode(mol, nm, mode_idx, Name, Value, ...)
%
%   Input:
%       mol       - struct returned by G16_structure (geometry)
%       nm        - struct returned by G16_nmodes (displacement vectors)
%       mode_idx  - mode index (1-based, indexing into nm.freq)
%
%   Optional parameters:
%       'Scale'       - arrow length scale (default: 1.5)
%       'ArrowColor'  - arrow colour (default: [1 0.4 0.1])
%       'AtomScale'   - CPK sphere scale factor (default: 0.35)
%       'BondTol'     - bond detection tolerance (default: 1.30)
%       'ShowLabels'  - show atom index labels (default: false)
%       'FlipSign'    - invert all arrow directions (default: false)
%                       Normal mode eigenvectors have an arbitrary overall
%                       sign (a mode and its 180°-phase-shifted twin are
%                       physically identical). If the arrows point opposite
%                       to GaussView's rendering for this mode, set this
%                       to true to flip them.
%
%   Example:
%       mol = G16_structure('V_E00t.out');
%       nm  = G16_nmodes('V_E00t.out');
%       G16_draw_mode(mol, nm, 91)                       % CC stretch a 1582 cm-1
%       G16_draw_mode(mol, nm, 91, 'FlipSign', true)      % reversed arrows

% -------------------------------------------------------------------------
p = inputParser;
addRequired(p, 'mol');
addRequired(p, 'nm');
addRequired(p, 'mode_idx', @isnumeric);
addParameter(p, 'Scale',      1.5,            @isnumeric);
addParameter(p, 'ArrowColor', [1.0 0.4 0.1],  @isnumeric);
addParameter(p, 'AtomScale',  0.35,           @isnumeric);
addParameter(p, 'BondTol',    1.30,           @isnumeric);
addParameter(p, 'ShowLabels', false,          @islogical);
addParameter(p, 'FlipSign',   false,          @islogical);
parse(p, mol, nm, mode_idx, varargin{:});

scale       = p.Results.Scale;
arrow_color = p.Results.ArrowColor;
atom_scale  = p.Results.AtomScale;
bond_tol    = p.Results.BondTol;
show_labels = p.Results.ShowLabels;
flip_sign   = p.Results.FlipSign;

% Validation
if mode_idx < 1 || mode_idx > nm.Nmodes
    error('G16_draw_mode: mode index %d is out of range [1, %d]', mode_idx, nm.Nmodes);
end
if mol.Natoms ~= nm.Natoms
    error('G16_draw_mode: mol.Natoms (%d) does not match nm.Natoms (%d)', mol.Natoms, nm.Natoms);
end

% -------------------------------------------------------------------------
% Render the CPK molecular structure
% -------------------------------------------------------------------------
freq_str = sprintf('Mode %d  —  %.1f cm^{-1}', mode_idx, nm.freq(mode_idx));
if nm.has_Raman
    freq_str = sprintf('%s   IR=%.1f   Raman=%.1f', ...
        freq_str, nm.IR(mode_idx), nm.Raman(mode_idx));
end

fig = figure('Color', 'white', 'Name', ...
    sprintf('Mode %d — %.1f cm-1', mode_idx, nm.freq(mode_idx)), ...
    'NumberTitle', 'off');
ax = axes('Parent', fig);

G16_draw_molecule(mol, 'Ax', ax, 'AtomScale', atom_scale, ...
    'BondTol', bond_tol, 'ShowLabels', show_labels, ...
    'ShowLegend', false, 'Title', freq_str);

% -------------------------------------------------------------------------
% Extract displacement vectors for the selected mode
% -------------------------------------------------------------------------
U = squeeze(nm.disp(:, :, mode_idx));   % [Natoms x 3]
if flip_sign
    U = -U;
end

% Normalise so that the largest per-atom displacement = 1 (scale-independent)
norms_i = sqrt(sum(U.^2, 2));     % per-atom Euclidean norm [Natoms x 1]
max_norm = max(norms_i);
if max_norm > 0
    U_scaled = U / max_norm * scale;
else
    warning('G16_draw_mode: zero displacement vectors for mode %d', mode_idx);
    return
end

% Draw arrows only for atoms whose displacement exceeds 5% of the maximum
thresh = 0.05;
hold(ax, 'on');

for i = 1:mol.Natoms
    if norms_i(i) / max_norm < thresh
        continue
    end
    x0 = mol.xyz(i,1);
    y0 = mol.xyz(i,2);
    z0 = mol.xyz(i,3);
    dx = U_scaled(i,1);
    dy = U_scaled(i,2);
    dz = U_scaled(i,3);

    draw_arrow3(ax, x0, y0, z0, dx, dy, dz, arrow_color);
end

% -------------------------------------------------------------------------
% Add title annotation
% -------------------------------------------------------------------------
[~, fname] = fileparts(mol.filename);
title(ax, {strrep(fname, '_', '\_'), freq_str}, ...
    'Interpreter', 'tex', 'FontSize', 10);

rotate3d(ax, 'on');

end  % G16_draw_mode


% =========================================================================
%  3D arrow with cone tip
% =========================================================================
function draw_arrow3(ax, x0, y0, z0, dx, dy, dz, color)
% Arrow shaft (line segment)
line(ax, [x0, x0+dx], [y0, y0+dy], [z0, z0+dz], ...
     'Color', color, 'LineWidth', 2.0, 'HandleVisibility', 'off');

% Cone tip aligned with the arrow direction
len  = sqrt(dx^2 + dy^2 + dz^2);
if len < 1e-6, return; end

tip_frac  = 0.25;   % cone occupies the last 25% of the total arrow length
tip_r     = 0.07;   % cone base radius in Angstroms

tip_len = len * tip_frac;
ux = dx/len; uy = dy/len; uz = dz/len;

% Start point of the cone (stepped back from the arrow tip)
cx = x0 + dx - ux*tip_len;
cy = y0 + dy - uy*tip_len;
cz = z0 + dz - uz*tip_len;

% Cone base: circle perpendicular to u
% Build a vector perpendicular to u (used as the cone base frame)
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

% Render the cone as a fan of triangular patches
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
