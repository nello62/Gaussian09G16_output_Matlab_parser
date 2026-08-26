function h = G_draw_density_surface(data, varargin)
% G_DRAW_DENSITY_SURFACE  Renders the total electron density (or, given a
%                          second structure, a density DIFFERENCE map; or,
%                          for an unrestricted calculation, the SPIN
%                          density alpha-minus-beta) as a real-space
%                          isosurface or 2D contour, evaluated from the
%                          raw .fchk basis-set data, overlaid on the CPK
%                          molecule -- a G_Utility equivalent of
%                          GaussView's "Total Density"/"Spin Density"
%                          surface types.
%
%   G_draw_density_surface(data)
%   G_draw_density_surface(data, 'Name', Value, ...)
%   G_draw_density_surface(data, 'CompareTo', data2, ...)
%   G_draw_density_surface(data, 'SpinDensity', true, ...)
%   h = G_draw_density_surface(...)
%
%   This is a standalone G_Utility companion function, independent of the
%   core toolbox (see G_DRAW_MO_SURFACE for the shared design rationale).
%   It shares its real-space Gaussian-type-orbital evaluation engine
%   (G_EVAL_MO_ON_GRID) with G_DRAW_MO_SURFACE: the total density is
%   ELECTRON-COUNT-EXACT, not a heuristic, since it is built from
%       rho(r) = sum_i n_i |psi_i(r)|^2
%   summed over all OCCUPIED molecular orbitals (n_i = 2 per orbital for
%   a closed-shell calculation), evaluated in a single pass through the
%   basis shells for every occupied orbital at once (see the "Method"
%   section below) rather than accumulating one orbital at a time.
%
%   Input:
%       data      - struct returned by G09_FCHK_READ or G16_FCHK_READ.
%                   Required fields: .filename, .mol, .Nbasis,
%                   .alpha_MO_coeff, .xyz_bohr, .Nalpha, .Nbeta. For the
%                   total density or 'CompareTo', only closed-shell
%                   (.Nalpha == .Nbeta) calculations are supported -- see
%                   the "Closed-shell only" note below. For
%                   'SpinDensity', an unrestricted (UHF/UKS) calculation
%                   is required instead (see 'SpinDensity' below) --
%                   .Nalpha == .Nbeta is NOT itself disqualifying (e.g. a
%                   broken-symmetry biradical singlet).
%
%   Optional parameters (Name-Value):
%       'CompareTo'    - a second data struct (same requirements as
%                        DATA). If given, renders the DIFFERENCE density
%                        rho(CompareTo) - rho(data) instead of the total
%                        density of data alone -- e.g. field-on minus
%                        field-off, to see where charge density
%                        increases/decreases. Both densities are
%                        evaluated on the same grid/plane, built from
%                        DATA's geometry (default: [], total density of
%                        data only). Mutually exclusive with
%                        'SpinDensity'
%       'SpinDensity'  - logical; if true, renders the SPIN density
%                        rho_alpha(r) - rho_beta(r) of DATA alone instead
%                        of the total density -- GaussView's "Spin
%                        Density" surface, showing where unpaired
%                        (alpha-excess, PosColor) vs.\ paired/beta-excess
%                        (NegColor) electron density is concentrated.
%                        Requires an unrestricted (UHF/UKS) calculation
%                        -- specifically, that the .fchk has a "Beta MO
%                        coefficients" section, re-read directly (like
%                        the basis-set data G_DRAW_MO_SURFACE already
%                        reads, no core-toolbox change); a restricted
%                        (RHF/RKS) job has identical alpha/beta orbitals
%                        by construction (zero spin density everywhere)
%                        and raises a clear error instead. Each
%                        spin-channel density uses occupation 1 per
%                        orbital (not 2, unlike the total/difference
%                        density case), since each alpha or beta
%                        spin-orbital holds exactly one electron.
%                        Mutually exclusive with 'CompareTo' (default:
%                        false)
%       'Mode'         - 'surface' (default): 3D isosurface(s). 'contour':
%                        a 2D filled-contour map on a single plane (see
%                        'Plane') -- as in G_DRAW_MO_SURFACE
%       'Plane'        - (only used when 'Mode' is 'contour') 'auto'
%                        (default) | 'xy' | 'xz' | 'yz' -- as in
%                        G_DRAW_MO_SURFACE
%       'PlaneOffset'  - (only used when 'Mode' is 'contour') shifts the
%                        evaluation plane along its own normal, in
%                        Angstrom (default: 0) -- as in G_DRAW_MO_SURFACE
%       'IsoValue'     - isosurface/contour-scale level, in electrons/Bohr^3.
%                        Default: 0.001 -- the classic Bader isodensity
%                        value conventionally used to define a molecule's
%                        outer envelope ("size"); appropriate for the
%                        TOTAL density (no 'CompareTo'). For a difference
%                        density, the relevant scale is usually smaller by
%                        an order of magnitude or more -- if no isosurface
%                        is found, a warning suggests reducing 'IsoValue'
%       'GridSpacing'  - real-space grid spacing, in Bohr (default: 0.15)
%       'Padding'      - grid/plane padding around the molecule's
%                        bounding box, in Bohr (default: 4.0); in
%                        'Mode','surface', auto-expanded up to 3 times if
%                        clipped, exactly as in G_DRAW_MO_SURFACE
%       'SaveCube'     - ('Mode','surface' only) filename; if given,
%                        writes the full volumetric grid actually used to
%                        build the isosurface (total density; density
%                        DIFFERENCE with 'CompareTo'; or SPIN density with
%                        'SpinDensity') to a Gaussian-format .cube file,
%                        at NO extra evaluation cost (it is the same grid
%                        already computed for the isosurface). Not
%                        available in 'Mode','contour' (default: '', no
%                        file written)
%       'PosColor'     - colour for the (only, without 'CompareTo'/
%                        'SpinDensity'; positive-difference or
%                        alpha-excess, with either) density lobe
%                        (default: [0.10 0.40 0.85], blue)
%       'NegColor'     - negative-difference (or beta-excess, with
%                        'SpinDensity') lobe colour, only used with
%                        'CompareTo'/'SpinDensity' (default:
%                        [0.85 0.15 0.10], red)
%       'FaceAlpha'    - lobe transparency, 0-1. ('Mode','surface' only.)
%                        Its own default (0.55) is used as-is with
%                        'SurfaceStyle','transparent'; overridden to 1.0
%                        with 'SurfaceStyle','solid' UNLESS given
%                        explicitly (which always wins); has no effect
%                        with 'SurfaceStyle','grid' (default: 0.55)
%       'SurfaceStyle' - ('Mode','surface' only) 'solid' | 'transparent'
%                        (default) | 'grid' -- as in GaussView's own
%                        surface-style choice. 'grid' draws the
%                        isosurface as a wireframe (mesh edges only, no
%                        face fill) coloured by 'PosColor'/'NegColor'
%                        instead of a filled patch -- useful for seeing
%                        through/inside a surface, or comparing overlaid
%                        surfaces, without stacking translucency. Since
%                        neither this function nor G_DRAW_MO_SURFACE
%                        decimates the 'Mode','surface' mesh (unlike
%                        G_DRAW_ESP_SURFACE's 'MaxVertices'), a fine
%                        'GridSpacing' can give a very dense wireframe --
%                        coarsen 'GridSpacing' for a readable grid render
%       'ShowMolecule' - overlay the molecule (3D CPK in 'surface' mode,
%                        a lightweight 2D sketch in 'contour' mode) --
%                        as in G_DRAW_MO_SURFACE (default: true)
%       'ShowLabels'   - as in G_DRAW_MO_SURFACE (default: false)
%       'AtomScale'    - as in G_DRAW_MO_SURFACE (default: 0.35)
%       'Title'        - figure title (default: auto)
%       'Ax'           - existing axes handle (default: new figure)
%       'Verbose'      - print grid size and timing info (default: true)
%
%   Output:
%       h - struct. In 'Mode','surface': .pos (patch handle for the
%           density, or the positive-difference/alpha-excess lobe with
%           'CompareTo'/'SpinDensity'; [] if not found) and .neg ([]
%           without 'CompareTo'/'SpinDensity'; the negative-difference/
%           beta-excess lobe patch handle otherwise). In 'Mode','contour':
%           .contour and .zero, as in G_DRAW_MO_SURFACE (.zero is always
%           [] without 'CompareTo'/'SpinDensity', since a total density
%           has no sign change to trace a nodal line through).
%
%   Method: reuses G_EVAL_MO_ON_GRID's real-space Gaussian-type-orbital
%   evaluator (see G_DRAW_MO_SURFACE for the full normalization/pure-
%   harmonic derivation), called ONCE with the full [Nbasis x Nocc]
%   occupied-MO-coefficient matrix (columns 1:Nalpha of the reshaped
%   alpha_MO_coeff, or 1:Nbeta of the reshaped beta MO coefficients for
%   the beta channel of 'SpinDensity'), returning an [Ngrid x Nocc]
%   matrix of every occupied orbital's amplitude at every grid point in
%   one pass through the basis shells (the per-shell radial/angular
%   evaluation, the expensive part, is shared across all occupied
%   orbitals rather than repeated once per orbital). The total/difference
%   density is rho(r) = 2 * sum_i psi_i(r)^2 (double occupancy); the spin
%   density is rho_alpha(r) - rho_beta(r), each with occupancy 1 (see
%   G_EVAL_DENSITY_ON_GRID's 'occ_factor' argument).
%
%   S, P, SP, D, F, and G shells (pure or Cartesian) are supported, with
%   the same validated normalization/ordering as G_DRAW_MO_SURFACE; a
%   basis set containing H (or higher) shells raises a clear error.
%
%   Closed-shell only (total density / 'CompareTo'): this function
%   requires data.Nalpha == data.Nbeta (and the same for 'CompareTo', if
%   given) and raises a clear error otherwise, rather than silently
%   building a wrong electron count for an open-shell (UHF/UKS)
%   calculation. 'SpinDensity' is the dedicated open-shell path instead
%   (see above) -- it does NOT gate on Nalpha==Nbeta (which is not itself
%   a restricted/unrestricted signal), gating instead on whether a "Beta
%   MO coefficients" .fchk section is actually present.
%
%   Example:
%       data = G16_fchk_read('4-NTP.fchk');
%       G_draw_density_surface(data);                          % total density
%       G_draw_density_surface(data, 'IsoValue', 0.002);
%
%       radical = G16_fchk_read('radical_uks.fchk');
%       G_draw_density_surface(radical, 'SpinDensity', true, 'IsoValue', 0.002);
%
%       data1 = G16_fchk_read('nofield.fchk');
%       data2 = G16_fchk_read('field_x025.fchk');
%       G_draw_density_surface(data1, 'CompareTo', data2, 'IsoValue', 0.0005);
%
%       G_draw_density_surface(data, 'SurfaceStyle', 'grid', 'GridSpacing', 0.3);
%
%   See also G_DRAW_MO_SURFACE, G09_FCHK_READ, G16_FCHK_READ.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

% -------------------------------------------------------------------------
% Parse input arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'data');
addParameter(p, 'CompareTo',    [],                @(s) isempty(s) || isstruct(s));
addParameter(p, 'SpinDensity',  false,             @islogical);
addParameter(p, 'Mode',         'surface',         @(s) any(strcmpi(s, {'surface','contour'})));
addParameter(p, 'Plane',        'auto',            @(s) any(strcmpi(s, {'auto','xy','xz','yz'})));
addParameter(p, 'PlaneOffset',  0,                 @isnumeric);
addParameter(p, 'IsoValue',     0.001,             @isnumeric);
addParameter(p, 'GridSpacing',  0.15,              @isnumeric);
addParameter(p, 'Padding',      4.0,               @isnumeric);
addParameter(p, 'SaveCube',     '',                @ischar);
addParameter(p, 'PosColor',     [0.10 0.40 0.85],  @isnumeric);
addParameter(p, 'NegColor',     [0.85 0.15 0.10],  @isnumeric);
addParameter(p, 'FaceAlpha',    0.55,              @isnumeric);
addParameter(p, 'SurfaceStyle', 'transparent',     @(s) any(strcmpi(s, {'solid','transparent','grid'})));
addParameter(p, 'ShowMolecule', true,              @islogical);
addParameter(p, 'ShowLabels',   false,             @islogical);
addParameter(p, 'AtomScale',    0.35,              @isnumeric);
addParameter(p, 'Title',        '',                @ischar);
addParameter(p, 'Ax',           [],                @ishandle);
addParameter(p, 'Verbose',      true,              @islogical);
parse(p, data, varargin{:});

data2        = p.Results.CompareTo;
compare_mode = ~isempty(data2);
spin_density = p.Results.SpinDensity;
signed_mode  = compare_mode || spin_density;
mode         = lower(p.Results.Mode);
plane_choice = lower(p.Results.Plane);
plane_offset = p.Results.PlaneOffset;
isoval       = abs(p.Results.IsoValue);
spacing      = p.Results.GridSpacing;
pad          = p.Results.Padding;
save_cube    = p.Results.SaveCube;
pos_color    = p.Results.PosColor;
neg_color    = p.Results.NegColor;
face_alpha   = p.Results.FaceAlpha;
surface_style = lower(p.Results.SurfaceStyle);
show_mol     = p.Results.ShowMolecule;
show_labels  = p.Results.ShowLabels;
atom_scale   = p.Results.AtomScale;
fig_title    = p.Results.Title;
ax           = p.Results.Ax;
verbose      = p.Results.Verbose;
padding_explicit = ~ismember('Padding', p.UsingDefaults);
face_alpha_explicit = ~ismember('FaceAlpha', p.UsingDefaults);

% 'SurfaceStyle' only sets FaceAlpha's DEFAULT (GaussView's "solid" =
% opaque, "transparent" = the existing 0.55 default); an explicit
% 'FaceAlpha' always wins, exactly as 'ShowMolecule' defers to an
% explicit 'FaceAlpha' in G_DRAW_ESP_SURFACE. 'grid' ignores FaceAlpha
% entirely (no face is drawn, only edges).
if ~face_alpha_explicit && strcmp(surface_style, 'solid')
    face_alpha = 1.0;
end
if face_alpha_explicit && strcmp(surface_style, 'grid') && verbose
    fprintf('G_draw_density_surface: ''FaceAlpha'' has no effect with ''SurfaceStyle'',''grid'' (no face is drawn, only mesh edges).\n');
end

if ~isempty(save_cube) && strcmp(mode, 'contour')
    error('G_draw_density_surface: ''SaveCube'' needs the full 3D grid built in ''Mode'',''surface'' -- not available in ''Mode'',''contour'' (a single 2D plane).');
end
if spin_density && compare_mode
    error('G_draw_density_surface: ''SpinDensity'' and ''CompareTo'' are mutually exclusive (spin density is alpha-minus-beta of DATA alone, not a comparison between two structs).');
end

% -------------------------------------------------------------------------
% Validate input struct(s) and build the occupied-MO coefficient matrix/
% matrices: closed-shell total/difference density needs only the alpha
% columns (doubly occupied); spin density needs the alpha AND beta
% columns separately (singly occupied each), from an unrestricted (UHF/
% UKS) calculation.
% -------------------------------------------------------------------------
if spin_density
    [aobasis1, occ_alpha, occ_beta] = validate_and_prepare_spin(data, 'data');
else
    [aobasis1, occ1] = validate_and_prepare(data, 'data');
    if compare_mode
        [aobasis2, occ2] = validate_and_prepare(data2, 'CompareTo');
    end
end

a0 = 0.529177210544;   % Bohr -> Angstrom (CODATA), consistent with the toolbox
xyz_bohr = data.xyz_bohr;

% -------------------------------------------------------------------------
% Set up figure/axes (shared by both modes)
% -------------------------------------------------------------------------
if isempty(ax)
    fig = figure('Color', 'white', 'Name', 'Density surface', 'NumberTitle', 'off');
    ax  = axes('Parent', fig);
end
hold(ax, 'on');

h.pos = [];
h.neg = [];
h.contour = [];
h.zero    = [];

if strcmp(mode, 'surface')
    % ---------------------------------------------------------------------
    % Build real-space grid and extract the isosurface(s), auto-expanding
    % Padding and retrying if clipped -- exactly as in G_draw_mo_surface.
    % ---------------------------------------------------------------------
    max_grid     = 30e6;
    max_attempts = 4;
    pad_growth   = 1.6;

    attempt = 0;
    clipped = false;
    while true
        attempt = attempt + 1;

        lo = min(xyz_bohr, [], 1) - pad;
        hi = max(xyz_bohr, [], 1) + pad;
        gx = lo(1):spacing:hi(1);
        gy = lo(2):spacing:hi(2);
        gz = lo(3):spacing:hi(3);
        Ngrid = numel(gx) * numel(gy) * numel(gz);

        if Ngrid > max_grid
            if attempt == 1
                error(['G_draw_density_surface: requested grid has %d points (%d x %d x %d), ' ...
                       'too large. Increase ''GridSpacing'' or decrease ''Padding''.'], ...
                      Ngrid, numel(gx), numel(gy), numel(gz));
            end
            warning('G_draw_density_surface: stopped auto-expanding ''Padding'' at %.1f Bohr to stay under the %d-point grid limit; the surface may still be clipped. Increase ''GridSpacing'' to allow more room.', pad, max_grid);
            break
        end

        if verbose
            fprintf('G_draw_density_surface: grid %d x %d x %d = %d points (spacing=%.2f Bohr, padding=%.1f Bohr)\n', ...
                numel(gx), numel(gy), numel(gz), Ngrid, spacing, pad);
        end

        [X, Y, Z] = meshgrid(gx, gy, gz);
        grid_pts_bohr = [X(:), Y(:), Z(:)];

        tic;
        if spin_density
            V = reshape(g_eval_density_on_grid(aobasis1, occ_alpha, grid_pts_bohr, 1), size(X)) - ...
                reshape(g_eval_density_on_grid(aobasis1, occ_beta,  grid_pts_bohr, 1), size(X));
        else
            V = reshape(g_eval_density_on_grid(aobasis1, occ1, grid_pts_bohr), size(X));
            if compare_mode
                V = reshape(g_eval_density_on_grid(aobasis2, occ2, grid_pts_bohr), size(X)) - V;
            end
        end
        if verbose
            fprintf('  Density evaluation: %.1f s\n', toc);
        end

        Xa = X * a0;
        Ya = Y * a0;
        Za = Z * a0;

        fv_pos = isosurface(Xa, Ya, Za, V, isoval);
        if signed_mode
            fv_neg = isosurface(Xa, Ya, Za, V, -isoval);
        else
            fv_neg = struct('vertices', [], 'faces', []);
        end

        clipped = g_is_clipped(fv_pos, Xa, Ya, Za, 1.5*spacing*a0) || ...
                  g_is_clipped(fv_neg, Xa, Ya, Za, 1.5*spacing*a0);

        if ~clipped || padding_explicit || attempt >= max_attempts
            break
        end

        if verbose
            fprintf('  Isosurface touches the grid boundary -- expanding Padding %.1f -> %.1f Bohr and retrying.\n', pad, pad*pad_growth);
        end
        pad = pad * pad_growth;
    end

    if clipped
        if padding_explicit
            warning('G_draw_density_surface: the isosurface appears to be clipped by the grid boundary (Padding=%.1f Bohr). Try a larger ''Padding''.', pad);
        else
            warning('G_draw_density_surface: the isosurface still appears clipped after auto-expanding ''Padding'' to %.1f Bohr (%d attempts). Try an even larger ''Padding'' explicitly, or a larger ''GridSpacing'' to afford more room within the grid-size limit.', pad, attempt);
        end
    end

    if ~isempty(save_cube)
        if spin_density
            title_line = 'Spin density (alpha - beta) -- G_Utility export';
        elseif compare_mode
            title_line = 'Density difference (CompareTo - data) -- G_Utility export';
        else
            title_line = 'Total electron density -- G_Utility export';
        end
        comment_line = sprintf('rho(r), electrons/Bohr^3, grid spacing=%.3f Bohr', spacing);
        g_write_cube_file(save_cube, title_line, comment_line, data.mol.Z, xyz_bohr, gx, gy, gz, V);
        if verbose
            fprintf('  Cube file written: %s\n', save_cube);
        end
    end

    % ---------------------------------------------------------------------
    % Draw isosurface(s)
    % ---------------------------------------------------------------------
    if ~isempty(fv_pos.vertices)
        pos_style_args = surface_style_args(surface_style, pos_color, face_alpha);
        h.pos = patch(ax, fv_pos, pos_style_args{:}, ...
            'FaceLighting', 'gouraud', 'DisplayName', '+');
        isonormals(Xa, Ya, Za, V, h.pos);
    end

    if signed_mode && ~isempty(fv_neg.vertices)
        neg_style_args = surface_style_args(surface_style, neg_color, face_alpha);
        h.neg = patch(ax, fv_neg, neg_style_args{:}, ...
            'FaceLighting', 'gouraud', 'DisplayName', '-');
        isonormals(Xa, Ya, Za, V, h.neg);
    end

    if isempty(fv_pos.vertices) && (~signed_mode || isempty(fv_neg.vertices))
        warning('G_draw_density_surface: no isosurface found at IsoValue=%.4g -- try a smaller ''IsoValue'' or larger ''Padding''.', isoval);
    end

    % ---------------------------------------------------------------------
    % Overlay molecule
    % ---------------------------------------------------------------------
    if show_mol
        if exist('G16_draw_molecule', 'file') == 2
            draw_mol_fcn = @G16_draw_molecule;
        elseif exist('G09_draw_molecule', 'file') == 2
            draw_mol_fcn = @G09_draw_molecule;
        else
            error('G_draw_density_surface: ''ShowMolecule'' is true but neither G16_draw_molecule nor G09_draw_molecule is on the MATLAB path.');
        end
        draw_mol_fcn(data.mol, 'Ax', ax, 'ShowLabels', show_labels, ...
            'AtomScale', atom_scale, 'ShowLegend', true, 'Title', '');
    else
        axis(ax, 'equal');
        axis(ax, 'off');
        view(ax, 3);
        lighting(ax, 'gouraud');
        material(ax, 'dull');
        camlight(ax, 'headlight');
        camlight(ax, 45, 30);
        axis(ax, 'tight');
        camproj(ax, 'perspective');
    end

else   % mode == 'contour'
    % ---------------------------------------------------------------------
    % Build the cutting plane and evaluate the density on a 2D grid within
    % it, exactly as in G_draw_mo_surface's 'contour' mode.
    % ---------------------------------------------------------------------
    [center, u, v, n] = g_plane_basis(xyz_bohr, plane_choice);
    offset_bohr = plane_offset / a0;

    atom_uv = (xyz_bohr - center) * [u, v];   % [Nat x 2], Bohr
    lo_uv = min(atom_uv, [], 1) - pad;
    hi_uv = max(atom_uv, [], 1) + pad;
    gs = lo_uv(1):spacing:hi_uv(1);
    gt = lo_uv(2):spacing:hi_uv(2);
    Ngrid = numel(gs) * numel(gt);

    max_grid = 4e6;
    if Ngrid > max_grid
        error(['G_draw_density_surface: requested contour grid has %d points (%d x %d), ' ...
               'too large. Increase ''GridSpacing'' or decrease ''Padding''.'], ...
              Ngrid, numel(gs), numel(gt));
    end

    if verbose
        fprintf('G_draw_density_surface: contour grid %d x %d = %d points (spacing=%.2f Bohr, padding=%.1f Bohr, plane=%s, offset=%.2f %c)\n', ...
            numel(gs), numel(gt), Ngrid, spacing, pad, plane_choice, plane_offset, 197);
    end

    [S, T] = meshgrid(gs, gt);
    grid_pts_bohr = (center + offset_bohr*n') + S(:)*u' + T(:)*v';

    tic;
    if spin_density
        V2 = reshape(g_eval_density_on_grid(aobasis1, occ_alpha, grid_pts_bohr, 1), size(S)) - ...
             reshape(g_eval_density_on_grid(aobasis1, occ_beta,  grid_pts_bohr, 1), size(S));
    else
        V2 = reshape(g_eval_density_on_grid(aobasis1, occ1, grid_pts_bohr), size(S));
        if compare_mode
            V2 = reshape(g_eval_density_on_grid(aobasis2, occ2, grid_pts_bohr), size(S)) - V2;
        end
    end
    if verbose
        fprintf('  Density evaluation: %.1f s\n', toc);
    end

    Sa = S * a0;
    Ta = T * a0;

    % The electron density has a very sharp cusp-like peak at each
    % nucleus (visible even on a plane that only grazes an atom), one to
    % several orders of magnitude larger than the chemically interesting
    % bonding/valence-region values. Scaling the color axis to the true
    % max would saturate the whole map to near-white everywhere except a
    % few nuclear pinpoints, hiding exactly the region of interest -- so
    % the color scale is capped at the 98th percentile of |V2| instead
    % (a standard fix for this well-known density-visualization dynamic-
    % range problem); values above the cap simply saturate to the
    % top-of-scale colour, via CONTOURF's normal out-of-range behaviour.
    if signed_mode
        % Diverging colormap (NegColor -> white -> PosColor), zero pinned
        % to white regardless of the actual data range.
        n_cmap = 256;
        cmap = [interp1([0 1], [neg_color; 1 1 1], linspace(0,1,n_cmap/2)); ...
                interp1([0 1], [1 1 1; pos_color], linspace(0,1,n_cmap/2))];
        colormap(ax, cmap);
        vmax = g_pctile_local(abs(V2(:)), 98);
        if vmax == 0, vmax = 1; end
        clow = -vmax; chigh = vmax;
    else
        % Sequential colormap (white -> PosColor); density is >= 0 everywhere.
        n_cmap = 256;
        cmap = interp1([0 1], [1 1 1; pos_color], linspace(0,1,n_cmap));
        colormap(ax, cmap);
        vmax = g_pctile_local(V2(:), 98);
        if vmax <= 0, vmax = 1; end
        clow = 0; chigh = vmax;
    end

    n_levels = 21;
    [~, h.contour] = contourf(ax, Sa, Ta, V2, linspace(clow, chigh, n_levels), 'LineStyle', 'none');
    caxis(ax, [clow, chigh]);   %#ok<CAXIS> -- kept over CLIM for MATLAB R2021b compatibility

    if signed_mode && vmax > 0 && any(V2(:) > 0, 'all') && any(V2(:) < 0, 'all')
        [~, h.zero] = contour(ax, Sa, Ta, V2, [0 0], 'Color', [0.15 0.15 0.15], 'LineWidth', 1.3);
    end

    cb = colorbar(ax);
    if spin_density
        cb.Label.String = 'Spin density, alpha - beta (e/Bohr^3)';
    elseif compare_mode
        cb.Label.String = 'Density difference (e/Bohr^3)';
    else
        cb.Label.String = 'Electron density (e/Bohr^3)';
    end

    axis(ax, 'equal');
    box(ax, 'on');
    xlabel(ax, sprintf('in-plane u (%c)', 197), 'Interpreter', 'none');
    ylabel(ax, sprintf('in-plane v (%c)', 197), 'Interpreter', 'none');

    if show_mol
        g_draw_2d_atoms(ax, data.mol.symbols, atom_uv * a0, xyz_bohr * a0);
    end
end

% -------------------------------------------------------------------------
% Title
% -------------------------------------------------------------------------
if isempty(fig_title)
    if spin_density
        base_str = 'Spin density (\alpha-\beta)';
    elseif compare_mode
        base_str = 'Density difference';
    else
        base_str = 'Total density';
    end
    if strcmp(mode, 'surface')
        fig_title = sprintf('%s, iso=%s%.3g', base_str, ternary(signed_mode, '\pm', ''), isoval);
    else
        fig_title = sprintf('%s, plane=%s', base_str, plane_choice);
    end
end
title(ax, fig_title, 'Interpreter', 'tex', 'FontSize', 11);

hold(ax, 'off');

end % G_draw_density_surface


% =========================================================================
%  Local functions
% =========================================================================

function [aobasis, occ_coeff] = validate_and_prepare(data, argname)
%VALIDATE_AND_PREPARE  Validates one data struct (either the primary DATA
%   or a 'CompareTo' struct), re-reads its raw basis-set data, checks it
%   is closed-shell, and returns the occupied-MO coefficient matrix
%   [Nbasis x Nalpha] (columns 1:Nalpha of the reshaped alpha_MO_coeff).
    required = {'filename', 'mol', 'Nbasis', 'Nbasis_indep', 'alpha_MO_coeff', 'xyz_bohr', 'Nalpha', 'Nbeta'};
    for k = 1:numel(required)
        if ~isfield(data, required{k})
            error('G_draw_density_surface: %s must be the struct returned by G09_fchk_read/G16_fchk_read (missing "%s" field).', argname, required{k});
        end
    end
    if ~isfile(data.filename)
        error('G_draw_density_surface: %s.filename (%s) not found -- this function re-reads the raw basis-set sections directly from the .fchk file.', argname, data.filename);
    end
    if data.Nalpha ~= data.Nbeta
        error('G_draw_density_surface: %s is an open-shell calculation (Nalpha=%d, Nbeta=%d) -- only closed-shell (Nalpha==Nbeta) total-density evaluation is currently supported.', argname, data.Nalpha, data.Nbeta);
    end

    aobasis = g_read_aobasis_from_fchk(data.filename);
    if isempty(aobasis.shell_types)
        error('G_draw_density_surface: no basis-set sections found in %s (%s) -- is this a valid .fchk file?', data.filename, argname);
    end
    g_check_supported_shells(aobasis.shell_types);

    alpha_MO = reshape(data.alpha_MO_coeff, data.Nbasis, data.Nbasis_indep);
    occ_coeff = alpha_MO(:, 1:data.Nalpha);
end

function out = ternary(cond, a, b)
%TERNARY  Small inline conditional helper for the title string.
    if cond
        out = a;
    else
        out = b;
    end
end

function args = surface_style_args(style, color, face_alpha)
%SURFACE_STYLE_ARGS  PATCH Name-Value pairs for 'SurfaceStyle': 'solid'/
%   'transparent' render the usual filled, unlit-edge surface (only
%   FaceAlpha differs, resolved by the caller before this is invoked);
%   'grid' instead draws GaussView-style wireframe/mesh -- no face fill,
%   just the isosurface triangulation's own edges, coloured COLOR.
    switch style
        case 'grid'
            args = {'FaceColor', 'none', 'EdgeColor', color, 'FaceAlpha', 1};
        otherwise   % 'solid' / 'transparent'
            args = {'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', face_alpha};
    end
end

function [aobasis, occ_alpha, occ_beta] = validate_and_prepare_spin(data, argname)
%VALIDATE_AND_PREPARE_SPIN  Like VALIDATE_AND_PREPARE, but for
%   'SpinDensity' mode: does NOT require Nalpha==Nbeta (spin density is
%   meaningful for any unrestricted calculation, including an
%   Nalpha==Nbeta broken-symmetry biradical singlet), and additionally
%   reads the "Beta MO coefficients" .fchk section directly (re-read
%   raw, exactly as G_READ_AOBASIS_FROM_FCHK re-reads the basis-set
%   sections, no core-toolbox change needed) -- present only for a
%   genuinely unrestricted (UHF/UKS) calculation; a restricted (RHF/RKS)
%   .fchk never writes this section (alpha and beta MOs are identical by
%   construction), so its absence is used here as the actual open-shell
%   test, not Nalpha~=Nbeta.
    required = {'filename', 'mol', 'Nbasis', 'Nbasis_indep', 'alpha_MO_coeff', 'xyz_bohr', 'Nalpha', 'Nbeta'};
    for k = 1:numel(required)
        if ~isfield(data, required{k})
            error('G_draw_density_surface: %s must be the struct returned by G09_fchk_read/G16_fchk_read (missing "%s" field).', argname, required{k});
        end
    end
    if ~isfile(data.filename)
        error('G_draw_density_surface: %s.filename (%s) not found -- this function re-reads the raw basis-set sections directly from the .fchk file.', argname, data.filename);
    end

    aobasis = g_read_aobasis_from_fchk(data.filename);
    if isempty(aobasis.shell_types)
        error('G_draw_density_surface: no basis-set sections found in %s (%s) -- is this a valid .fchk file?', data.filename, argname);
    end
    g_check_supported_shells(aobasis.shell_types);

    alpha_MO = reshape(data.alpha_MO_coeff, data.Nbasis, data.Nbasis_indep);
    occ_alpha = alpha_MO(:, 1:data.Nalpha);

    beta_vec = read_beta_mo_coeff(data.filename, data.Nbasis, data.Nbasis_indep);
    if isempty(beta_vec)
        error('G_draw_density_surface: ''SpinDensity'' requires an unrestricted (UHF/UKS) calculation, but no "Beta MO coefficients" section was found in %s (%s) -- a restricted (RHF/RKS) job has identical alpha/beta orbitals and zero spin density everywhere.', data.filename, argname);
    end
    beta_MO = reshape(beta_vec, data.Nbasis, data.Nbasis_indep);
    occ_beta = beta_MO(:, 1:data.Nbeta);
end

function vals = read_beta_mo_coeff(filename, Nbasis, Nbasis_indep)
%READ_BETA_MO_COEFF  Reads the .fchk "Beta MO coefficients" section
%   directly (a minimal, single-purpose parser, in the same
%   self-contained spirit as G_READ_AOBASIS_FROM_FCHK/
%   G_draw_esp_surface's READ_NUCLEAR_CHARGES) -- present only for an
%   unrestricted (UHF/UKS) calculation. Returns [] (not an error) if the
%   section is absent, so the caller can distinguish "restricted
%   calculation" from a malformed file. VALS is the flat
%   [Nbasis*Nbasis_indep x 1] vector, same column-major convention as
%   data.alpha_MO_coeff (reshape(vals, Nbasis, Nbasis_indep) gives the
%   coefficient matrix, columns = MOs; Nbasis_indep may be less than
%   Nbasis if Gaussian dropped near-linear-dependent combinations).
    fid = fopen(filename, 'r');
    raw = fread(fid, '*char')';
    fclose(fid);
    lines = strsplit(raw, newline);
    N = numel(lines);

    header_re = '^(.{1,43}?)\s{1,3}(I|R|C)\s+N=\s*(\d+)';
    for k = 1:N
        ln = lines{k};
        if numel(ln) < 45, continue; end
        tok = regexp(ln, header_re, 'tokens', 'once');
        if isempty(tok), continue; end
        if ~strcmpi(strtrim(tok{1}), 'Beta MO coefficients'), continue; end
        nvals = str2double(tok{3});
        v = [];
        k2 = k + 1;
        while numel(v) < nvals && k2 <= N
            ln2 = strtrim(lines{k2});
            if isempty(ln2), k2 = k2+1; continue; end
            if numel(ln2) > 44 && ~isempty(regexp(ln2, header_re, 'once'))
                break
            end
            v = [v; sscanf(ln2, '%f')]; %#ok<AGROW>
            k2 = k2 + 1;
        end
        vals = v(1:min(end,nvals));
        vals = vals(:);
        expected = Nbasis * Nbasis_indep;
        if numel(vals) ~= expected
            error('G_draw_density_surface: "Beta MO coefficients" section in %s has %d entries, expected %d (Nbasis*Nbasis_indep).', filename, numel(vals), expected);
        end
        return
    end
    vals = [];   % section not found -- restricted calculation
end

