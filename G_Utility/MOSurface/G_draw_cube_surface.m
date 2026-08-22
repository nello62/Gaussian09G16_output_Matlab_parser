function h = G_draw_cube_surface(cubefile, varargin)
% G_DRAW_CUBE_SURFACE  Renders an isosurface directly from a saved
%                       Gaussian-format .cube file -- no .fchk, no
%                       basis-set re-evaluation. The load-and-replot
%                       companion to G_DRAW_DENSITY_SURFACE/
%                       G_DRAW_MO_SURFACE/G_DRAW_ESP_SURFACE's own
%                       'SaveCube' option (and works equally well on a
%                       real cubegen/GaussView-generated .cube).
%
%   G_draw_cube_surface(cubefile)
%   G_draw_cube_surface(cubefile, 'Name', Value, ...)
%   G_draw_cube_surface(cubefile, 'ColorBy', esp_cubefile, ...)
%   h = G_draw_cube_surface(...)
%
%   Input:
%       cubefile  - path to a .cube file (e.g. one written by any of the
%                   three G_DRAW_*_SURFACE functions' 'SaveCube' option,
%                   or by Gaussian's own cubegen).
%
%   Optional parameters (Name-Value):
%       'ColorBy'      - path to a SECOND .cube file. If given, the
%                         isosurface SHAPE still comes from CUBEFILE (a
%                         single positive-value envelope, as for a
%                         density cube), but each vertex is coloured by
%                         trilinearly interpolating the SECOND cube's
%                         field at that point -- reproducing GaussView's
%                         classic "ESP mapped on electron density"
%                         picture from two independently-saved cube
%                         files (a density cube as CUBEFILE, an ESP cube
%                         as 'ColorBy'), with no .fchk needed at all.
%                         The two cubes must share an IDENTICAL grid
%                         (same origin, spacing, and point counts -- use
%                         the same 'GridSpacing'/'Padding' when saving
%                         both from G_DRAW_DENSITY_SURFACE/
%                         G_DRAW_ESP_SURFACE's 'SaveCube', or matching
%                         'CubeSpacing'/'CubePadding' for the ESP one);
%                         a clear error is raised otherwise, rather than
%                         silently resampling (default: '', no coloring
%                         cube; CUBEFILE's own field is both shape and,
%                         if signed, colour -- see 'IsoValue' below)
%       'IsoValue'     - isosurface level, in the cube's own units.
%                         Default: guessed from CUBEFILE's title/comment
%                         text -- 0.001 for a density-like field (title
%                         contains "density"), 0.02 for an MO-amplitude
%                         cube (negative-Natoms header, GaussView's own
%                         MO default), matching G_DRAW_DENSITY_SURFACE/
%                         G_DRAW_MO_SURFACE's own defaults for a cube
%                         saved by this toolbox. For anything else
%                         (an ESP cube, or an unrecognized/external
%                         cube's title) there is no universally sensible
%                         default -- 'IsoValue' must then be given
%                         explicitly, or a clear error explains why
%       'SurfaceStyle' - 'solid' | 'transparent' (default) | 'grid' --
%                         as in the three G_DRAW_*_SURFACE functions.
%                         'grid' draws a wireframe; with 'ColorBy', its
%                         edges are coloured by interpolating the second
%                         cube's field (EdgeColor,'interp'), like
%                         G_DRAW_ESP_SURFACE; without 'ColorBy', by
%                         'PosColor'/'NegColor' like the others
%       'FaceAlpha'    - surface opacity, 0-1. Used as-is with
%                         'SurfaceStyle','transparent' (default: 0.55);
%                         overridden to 1.0 with 'solid' unless given
%                         explicitly (which always wins); unused with
%                         'grid'
%       'PosColor'     - colour for the positive lobe (no 'ColorBy'), or
%                         the positive end of the diverging colormap
%                         (with 'ColorBy') (default: [0.10 0.40 0.85], blue)
%       'NegColor'     - colour for the negative lobe/colormap end, only
%                         relevant for a SIGNED field (default:
%                         [0.85 0.15 0.10], red) -- see "Method" below
%                         for how signed-ness is decided
%       'MaxVertices'  - the isosurface mesh is decimated (MATLAB's
%                         built-in REDUCEPATCH) to at most this many
%                         vertices, mainly to bound 'ColorBy' interpolation
%                         cost and keep the render light (default: 5000)
%       'ShowMolecule' - overlay a CPK ball-and-stick model built
%                         directly from the atoms stored IN the cube
%                         file (no .fchk needed) (default: true)
%       'ShowLabels'   - as in G_DRAW_MO_SURFACE (default: false)
%       'AtomScale'    - as in G_DRAW_MO_SURFACE (default: 0.35)
%       'Title'        - figure title (default: auto, from the cube's
%                         own title line)
%       'Ax'           - existing axes handle (default: new figure)
%       'Verbose'      - print grid size and field-type detection info
%                         (default: true)
%
%   Output:
%       h - struct with fields .pos/.neg (patch handles for the
%           positive/negative lobe; .neg is [] for an unsigned field or
%           when 'ColorBy' is given, since that mode always renders a
%           single shape-defining lobe).
%
%   Method: unlike the three G_DRAW_*_SURFACE functions (which evaluate
%   Gaussian-type orbitals on a freshly built grid), this function reads
%   an already-computed scalar field straight from CUBEFILE via
%   G_READ_CUBE_FILE and extracts its isosurface(s) directly with
%   MATLAB's ISOSURFACE -- no basis-set data, no .fchk, no re-evaluation
%   of anything. A field is treated as SIGNED (drawing both a +IsoValue
%   and a -IsoValue lobe, coloured 'PosColor'/'NegColor') if it is an MO
%   cube (negative Natoms header) or if its values genuinely span both
%   signs beyond numerical noise; otherwise (e.g. a total density) only
%   the single +IsoValue lobe is drawn. 'ColorBy' overrides this: the
%   shape is always the single +IsoValue envelope of CUBEFILE, coloured
%   by the second cube's interpolated value instead of a flat colour.
%
%   Example:
%       % Re-plot a density cube saved earlier, no .fchk needed:
%       G_draw_cube_surface('density.cube');
%
%       % Re-plot a saved MO cube:
%       G_draw_cube_surface('homo.cube', 'SurfaceStyle', 'grid');
%
%       % Classic "ESP mapped on density" from two independently saved
%       % cubes (see G_DRAW_DENSITY_SURFACE/G_DRAW_ESP_SURFACE's
%       % 'SaveCube'), entirely without re-reading the .fchk:
%       G_draw_cube_surface('density.cube', 'ColorBy', 'esp.cube', ...
%           'IsoValue', 0.001);
%
%   See also G_DRAW_DENSITY_SURFACE, G_DRAW_MO_SURFACE, G_DRAW_ESP_SURFACE.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

% -------------------------------------------------------------------------
% Parse input arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'cubefile', @ischar);
addParameter(p, 'ColorBy',      '',                @ischar);
addParameter(p, 'IsoValue',     [],                @(x) isempty(x) || isnumeric(x));
addParameter(p, 'SurfaceStyle', 'transparent',     @(s) any(strcmpi(s, {'solid','transparent','grid'})));
addParameter(p, 'FaceAlpha',    0.55,              @isnumeric);
addParameter(p, 'PosColor',     [0.10 0.40 0.85],  @isnumeric);
addParameter(p, 'NegColor',     [0.85 0.15 0.10],  @isnumeric);
addParameter(p, 'MaxVertices',  5000,              @isnumeric);
addParameter(p, 'ShowMolecule', true,              @islogical);
addParameter(p, 'ShowLabels',   false,             @islogical);
addParameter(p, 'AtomScale',    0.35,              @isnumeric);
addParameter(p, 'Title',        '',                @ischar);
addParameter(p, 'Ax',           [],                @ishandle);
addParameter(p, 'Verbose',      true,              @islogical);
parse(p, cubefile, varargin{:});

color_by      = p.Results.ColorBy;
isoval_in     = p.Results.IsoValue;
surface_style = lower(p.Results.SurfaceStyle);
face_alpha    = p.Results.FaceAlpha;
face_alpha_explicit = ~ismember('FaceAlpha', p.UsingDefaults);
pos_color     = p.Results.PosColor;
neg_color     = p.Results.NegColor;
max_vertices  = p.Results.MaxVertices;
show_mol      = p.Results.ShowMolecule;
show_labels   = p.Results.ShowLabels;
atom_scale    = p.Results.AtomScale;
fig_title     = p.Results.Title;
ax            = p.Results.Ax;
verbose       = p.Results.Verbose;
color_mode    = ~isempty(color_by);

if ~face_alpha_explicit && strcmp(surface_style, 'solid')
    face_alpha = 1.0;
end
if face_alpha_explicit && strcmp(surface_style, 'grid') && verbose
    fprintf('G_draw_cube_surface: ''FaceAlpha'' has no effect with ''SurfaceStyle'',''grid'' (no face is drawn, only mesh edges).\n');
end

if ~isfile(cubefile)
    error('G_draw_cube_surface: cubefile (%s) not found.', cubefile);
end

% -------------------------------------------------------------------------
% Read the primary cube
% -------------------------------------------------------------------------
[atomic_numbers, xyz_bohr, gx, gy, gz, V, orbital_indices, title_line] = g_read_cube_file(cubefile);
is_mo_cube = ~isempty(orbital_indices);

if verbose
    fprintf('G_draw_cube_surface: %s -- grid %d x %d x %d, %d atoms%s\n', ...
        cubefile, numel(gx), numel(gy), numel(gz), numel(atomic_numbers), ...
        ternary(is_mo_cube, sprintf(' (MO cube, orbital %s)', mat2str(orbital_indices)), ''));
end

% -------------------------------------------------------------------------
% IsoValue: explicit value wins; otherwise guess from the title text
% written by this toolbox's own 'SaveCube' option (density-like -> 0.001,
% MO -> 0.02); anything else (an ESP cube, or an external/unrecognized
% cube) has no sensible universal default and must be given explicitly.
% -------------------------------------------------------------------------
if ~isempty(isoval_in)
    isoval = abs(isoval_in);
elseif is_mo_cube
    isoval = 0.02;
    if verbose, fprintf('  ''IsoValue'' not given -- guessed 0.02 (MO cube).\n'); end
elseif contains(lower(title_line), 'density')
    isoval = 0.001;
    if verbose, fprintf('  ''IsoValue'' not given -- guessed 0.001 (density-like title: "%s").\n', title_line); end
else
    error(['G_draw_cube_surface: ''IsoValue'' was not given, and the field type could not be ' ...
           'guessed from the cube''s title ("%s") -- there is no universal default for an ' ...
           'arbitrary/ESP-like field. Pass ''IsoValue'' explicitly.'], title_line);
end

% A lone ESP cube isosurfaced directly (no 'ColorBy') gives the surface
% of CONSTANT ESP VALUE -- a fundamentally different picture from
% "ESP mapped on the density envelope" (G_DRAW_ESP_SURFACE's own
% rendering, and the classic GaussView/chemist's mental image), which
% needs the DENSITY field for shape and the ESP field only for colour.
% No 'IsoValue' choice here reproduces that picture, because the
% shape-defining field itself differs -- flag this explicitly rather
% than let a "0.001" guess (that value's natural home is density,
% electrons/Bohr^3, not ESP, Hartree/e) silently produce a misleadingly
% "strange" surface with no explanation.
if ~color_mode && (contains(lower(title_line), 'electrostatic potential') || contains(lower(title_line), 'esp '))
    warning(['G_draw_cube_surface: %s appears to be an ESP field (title: "%s"), isosurfaced ' ...
             'here directly at IsoValue=%.4g -- this is a surface of CONSTANT ESP VALUE, not ' ...
             '"ESP mapped on the density envelope" (that needs the DENSITY cube for shape and ' ...
             'this ESP cube only for colour, via ''ColorBy''; see G_DRAW_CUBE_SURFACE''s help).'], ...
            cubefile, title_line, isoval);
end

% -------------------------------------------------------------------------
% Signed-ness: MO cube -> always signed; 'ColorBy' -> always a single
% shape-defining positive lobe (matching the classic density-envelope
% convention); otherwise auto-detect from the data itself.
% -------------------------------------------------------------------------
if color_mode
    signed_field = false;
elseif is_mo_cube
    signed_field = true;
else
    vtol = 1e-6 * max(abs(V(:)));
    signed_field = any(V(:) > vtol) && any(V(:) < -vtol);
end

a0 = 0.529177210544;   % Bohr -> Angstrom (CODATA), consistent with the toolbox
Xa = gx * a0; Ya = gy * a0; Za = gz * a0;
[X, Y, Z] = meshgrid(Xa, Ya, Za);

fv_pos = isosurface(X, Y, Z, V, isoval);
if signed_field
    fv_neg = isosurface(X, Y, Z, V, -isoval);
else
    fv_neg = struct('vertices', [], 'faces', []);
end

if isempty(fv_pos.vertices) && (~signed_field || isempty(fv_neg.vertices))
    error('G_draw_cube_surface: no isosurface found at IsoValue=%.4g -- the cube''s value range is [%.4g, %.4g]. Try a different ''IsoValue''.', ...
        isoval, min(V(:)), max(V(:)));
end

% -------------------------------------------------------------------------
% Decimate large meshes (mainly to bound 'ColorBy' interpolation cost)
% -------------------------------------------------------------------------
if size(fv_pos.vertices,1) > max_vertices
    n0 = size(fv_pos.vertices,1);
    fv_pos = reducepatch(fv_pos, max_vertices/n0);
    if verbose, fprintf('  Positive lobe decimated: %d -> %d vertices\n', n0, size(fv_pos.vertices,1)); end
end
if signed_field && size(fv_neg.vertices,1) > max_vertices
    n0 = size(fv_neg.vertices,1);
    fv_neg = reducepatch(fv_neg, max_vertices/n0);
    if verbose, fprintf('  Negative lobe decimated: %d -> %d vertices\n', n0, size(fv_neg.vertices,1)); end
end

% -------------------------------------------------------------------------
% 'ColorBy': read the second cube, require an identical grid, and
% trilinearly interpolate its field at every (decimated) vertex.
% -------------------------------------------------------------------------
if color_mode
    if ~isfile(color_by)
        error('G_draw_cube_surface: ''ColorBy'' file (%s) not found.', color_by);
    end
    [~, ~, gx2, gy2, gz2, V2] = g_read_cube_file(color_by);
    if numel(gx2)~=numel(gx) || numel(gy2)~=numel(gy) || numel(gz2)~=numel(gz) || ...
       max(abs(gx2(:)-gx(:)),[],'all') > 1e-4 || max(abs(gy2(:)-gy(:)),[],'all') > 1e-4 || max(abs(gz2(:)-gz(:)),[],'all') > 1e-4
        error(['G_draw_cube_surface: ''ColorBy'' (%s) does not share an identical grid with %s ' ...
               '(origin/spacing/point-counts must match exactly). Re-save both cubes with the ' ...
               'same grid settings (''GridSpacing''/''Padding'' or ''CubeSpacing''/''CubePadding'').'], ...
              color_by, cubefile);
    end
    verts_bohr = fv_pos.vertices / a0;
    cdata = interp3(X/a0, Y/a0, Z/a0, V2, verts_bohr(:,1), verts_bohr(:,2), verts_bohr(:,3));
    if verbose
        fprintf('  ''ColorBy'' field interpolated at %d vertices: range [%.4g, %.4g]\n', numel(cdata), min(cdata), max(cdata));
    end
end

% -------------------------------------------------------------------------
% Draw
% -------------------------------------------------------------------------
if isempty(ax)
    fig = figure('Color', 'white', 'Name', 'Cube surface', 'NumberTitle', 'off');
    ax  = axes('Parent', fig);
end
hold(ax, 'on');

h.pos = [];
h.neg = [];

if color_mode
    n_cmap = 256;
    cmap = [interp1([0 1], [neg_color; 1 1 1], linspace(0,1,n_cmap/2)); ...
            interp1([0 1], [1 1 1; pos_color], linspace(0,1,n_cmap/2))];
    colormap(ax, cmap);
    vmax = g_pctile_local(abs(cdata), 99);   % matches G_DRAW_ESP_SURFACE's own
                                              % colour-scale convention: the 99th
                                              % percentile, not the raw max, so a
                                              % single outlier vertex doesn't wash
                                              % out the whole colour scale
    if vmax == 0, vmax = 1; end

    if strcmp(surface_style, 'grid')
        cb_style_args = {'FaceColor', 'none', 'EdgeColor', 'interp', 'FaceAlpha', 1};
    else
        cb_style_args = {'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', face_alpha};
    end
    h.pos = patch(ax, 'Vertices', fv_pos.vertices, 'Faces', fv_pos.faces, ...
        'FaceVertexCData', cdata, cb_style_args{:}, 'FaceLighting', 'gouraud');
    caxis(ax, [-vmax, vmax]);   %#ok<CAXIS> -- kept over CLIM for MATLAB R2021b compatibility
    cb = colorbar(ax);
    cb.Label.String = 'ColorBy field value';
else
    if ~isempty(fv_pos.vertices)
        pos_args = surface_style_args(surface_style, pos_color, face_alpha);
        h.pos = patch(ax, fv_pos, pos_args{:}, 'FaceLighting', 'gouraud', 'DisplayName', '+');
        isonormals(X, Y, Z, V, h.pos);
    end
    if signed_field && ~isempty(fv_neg.vertices)
        neg_args = surface_style_args(surface_style, neg_color, face_alpha);
        h.neg = patch(ax, fv_neg, neg_args{:}, 'FaceLighting', 'gouraud', 'DisplayName', '-');
        isonormals(X, Y, Z, V, h.neg);
    end
end

axis(ax, 'equal');
axis(ax, 'off');
view(ax, 3);
lighting(ax, 'gouraud');
material(ax, 'dull');
camlight(ax, 'headlight');
camlight(ax, 45, 30);

if show_mol
    mol.Natoms  = numel(atomic_numbers);
    mol.symbols = arrayfun(@atomic_symbol, atomic_numbers, 'UniformOutput', false);
    mol.xyz     = xyz_bohr * a0;
    mol.filename = cubefile;
    if exist('G16_draw_molecule', 'file') == 2
        draw_mol_fcn = @G16_draw_molecule;
    elseif exist('G09_draw_molecule', 'file') == 2
        draw_mol_fcn = @G09_draw_molecule;
    else
        error('G_draw_cube_surface: ''ShowMolecule'' is true but neither G16_draw_molecule nor G09_draw_molecule is on the MATLAB path.');
    end
    draw_mol_fcn(mol, 'Ax', ax, 'ShowLabels', show_labels, ...
        'AtomScale', atom_scale, 'ShowLegend', true, 'Title', '');
end

axis(ax, 'tight');
camproj(ax, 'perspective');

if isempty(fig_title)
    fig_title = sprintf('%s (iso=%.3g)', title_line, isoval);
end
title(ax, fig_title, 'Interpreter', 'none', 'FontSize', 11);

hold(ax, 'off');

end % G_draw_cube_surface


% =========================================================================
%  Local functions
% =========================================================================

function out = ternary(cond, a, b)
%TERNARY  Small inline conditional helper.
    if cond
        out = a;
    else
        out = b;
    end
end

function args = surface_style_args(style, color, face_alpha)
%SURFACE_STYLE_ARGS  PATCH Name-Value pairs for 'SurfaceStyle' -- see
%   G_DRAW_DENSITY_SURFACE/G_DRAW_MO_SURFACE for the identical helper.
    switch style
        case 'grid'
            args = {'FaceColor', 'none', 'EdgeColor', color, 'FaceAlpha', 1};
        otherwise
            args = {'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', face_alpha};
    end
end

function sym = atomic_symbol(z)
%ATOMIC_SYMBOL  Atomic number -> element symbol, self-contained (same
%   table as G16_FCHK_READ's, duplicated here so this file never needs
%   the core toolbox -- consistent with every other G_Utility function).
    SYM = { ...
        'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', ...
        'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca', ...
        'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', ...
        'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr', ...
        'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', ...
        'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', ...
        'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', ...
        'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', ...
        'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', ...
        'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', ...
        'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', ...
        'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'};
    if z >= 1 && z <= numel(SYM)
        sym = SYM{z};
    else
        sym = sprintf('Z%d', z);
    end
end
