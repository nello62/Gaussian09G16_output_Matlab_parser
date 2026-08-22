function h = G_draw_mo_surface(data, mo_index, varargin)
% G_DRAW_MO_SURFACE  Renders a molecular-orbital real-space isosurface
%                     (positive/negative lobes), evaluated from the raw
%                     .fchk basis-set data, overlaid on the CPK molecule
%                     -- a G_Utility equivalent of GaussView's "Molecular
%                     Orbitals" surface type.
%
%   G_draw_mo_surface(data, mo_index)
%   G_draw_mo_surface(data, mo_index, Name, Value, ...)
%   h = G_draw_mo_surface(...)
%
%   This is a standalone G_Utility companion function, independent of the
%   core toolbox: it does NOT require (and does not use) any modification
%   to G09_FCHK_READ/G16_FCHK_READ. The handful of extra raw basis-set
%   .fchk sections needed to evaluate a molecular orbital in real space
%   (shell types, primitive exponents, contraction coefficients, ...) are
%   re-read directly from data.filename by this function's own minimal
%   parser (see local function READ_AOBASIS_FROM_FCHK below) -- so it
%   plugs into the existing, unmodified toolbox exactly like the other
%   G_Utility tools (e.g. G_STRUCTURAL_POLARIZABILITY) that consume a
%   G09/G16_fchk_read struct without requiring any new fields from it.
%
%   Input:
%       data      - struct returned by G09_FCHK_READ or G16_FCHK_READ.
%                   Required fields: .filename, .mol, .Nbasis,
%                   .alpha_MO_coeff, .xyz_bohr. Optional: .HOMO_idx,
%                   .alpha_orb_energies (used only for 'HOMO'/'LUMO'
%                   index resolution and the default title).
%       mo_index  - which MO to render, either:
%                     * a positive integer (1-based column of
%                       reshape(data.alpha_MO_coeff, data.Nbasis, data.Nbasis))
%                     * 'HOMO', 'LUMO', 'HOMO-n', 'LUMO+n' (n = integer),
%                       resolved against data.HOMO_idx
%
%   Optional parameters (Name-Value):
%       'Mode'         - 'surface' (default): 3D +/-isosurface, as above.
%                        'contour': a 2D filled-contour amplitude map on a
%                        single plane through the molecule (see 'Plane'),
%                        with the nodal (zero-amplitude) line drawn
%                        explicitly -- much cheaper than a 3D grid, and
%                        often clearer for reading off the nodal pattern
%                        of a planar/near-planar pi system than a 3D
%                        isosurface viewed from an angle.
%       'Plane'        - (only used when 'Mode' is 'contour') 'auto'
%                        (default): the best-fit plane through all atoms,
%                        found via PCA (the plane minimizing the atoms'
%                        summed squared out-of-plane distance) -- the
%                        natural choice for a planar/near-planar
%                        molecule. 'xy'/'xz'/'yz': a coordinate plane
%                        through the molecule's centroid.
%       'PlaneOffset'  - (only used when 'Mode' is 'contour') shifts the
%                        evaluation plane along its own normal by this
%                        many Angstrom (default: 0, i.e. the plane
%                        itself). Essential for a pi-symmetry orbital:
%                        its amplitude is antisymmetric about (i.e.
%                        exactly zero on) the molecular plane, so the
%                        default offset=0 shows nothing for a typical
%                        aromatic HOMO/LUMO -- use e.g. 0.5-1.0 to probe
%                        a parallel plane where the pi density peaks. The
%                        atom sketch is always drawn at the true,
%                        unshifted in-plane positions, regardless of this
%                        offset, for reference. A warning is printed if
%                        the resulting amplitude is essentially zero
%                        everywhere, suggesting this is the cause.
%       'IsoValue'     - wavefunction-amplitude isosurface level, in
%                        atomic units (default: 0.02, GaussView's own
%                        default). Both +IsoValue and -IsoValue lobes are
%                        drawn. ('Mode','surface' only.)
%       'GridSpacing'  - real-space grid spacing, in Bohr (default: 0.15)
%       'Padding'      - grid/plane padding around the molecule's
%                        bounding box, in Bohr (default: 4.0). In
%                        'Mode','surface', auto-expanded up to 3 times,
%                        x1.6 each time, if the isosurface touches the
%                        grid boundary -- see "Automatic padding" below
%                        (passing an explicit value disables this and a
%                        warning is issued instead if still clipped).
%                        Not auto-expanded in 'Mode','contour' (a
%                        cropped-looking 2D map is a normal, expected
%                        edge effect, not a hidden artifact the way a
%                        clipped 3D isosurface facet is).
%       'SaveCube'     - ('Mode','surface' only) filename; if given,
%                        writes the full volumetric grid actually used to
%                        build the isosurface (the MO amplitude, signed)
%                        to a Gaussian-format MO .cube file (negative
%                        Natoms + the orbital-index line, matching the
%                        format read/validated against real cubegen MO
%                        cubes), at NO extra evaluation cost. Not
%                        available in 'Mode','contour' (default: '', no
%                        file written)
%       'PosColor'     - positive-amplitude colour (default: [0.10 0.40 0.85], blue)
%       'NegColor'     - negative-amplitude colour (default: [0.85 0.15 0.10], red)
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
%                        through/inside a lobe, or comparing overlaid
%                        surfaces, without stacking translucency. This
%                        function does not decimate the 'Mode','surface'
%                        mesh (unlike G_DRAW_ESP_SURFACE's
%                        'MaxVertices'), so a fine 'GridSpacing' can give
%                        a very dense wireframe -- coarsen 'GridSpacing'
%                        for a readable grid render
%       'ShowMolecule' - overlay the CPK ball-and-stick model via
%                        G09/G16_DRAW_MOLECULE, whichever is found on the
%                        MATLAB path. In 'Mode','contour', overlays a
%                        lightweight 2D atom/bond sketch instead (no
%                        toolbox call). (default: true)
%       'ShowLabels'   - atom labels, forwarded to the molecule drawer in
%                        'Mode','surface' (default: false, usually
%                        cluttered with a surface); always shown as
%                        element symbols in 'Mode','contour'
%       'AtomScale'    - forwarded to the molecule drawer. ('Mode','surface' only.) (default: 0.35)
%       'Title'        - figure title (default: auto, includes MO index,
%                        HOMO/LUMO tag if applicable, and its energy)
%       'Ax'           - existing axes handle (default: new figure)
%       'Verbose'      - print grid size and timing info (default: true)
%
%   Output:
%       h - struct. In 'Mode','surface': .pos/.neg, the positive/negative
%           lobe patch handles ([] if that lobe was not found at
%           IsoValue). In 'Mode','contour': .contour (filled-contour
%           handle) and .zero (nodal-line handle, [] if the amplitude
%           does not change sign anywhere on the plane).
%
%   Method: each contracted Gaussian-type basis function is evaluated
%   directly on a real-space rectangular grid (primitive normalization
%   N(alpha,l,m,n) = (2*alpha/pi)^(3/4) * sqrt((4*alpha)^(l+m+n) /
%   ((2l-1)!!(2m-1)!!(2n-1)!!)); pure D/F/G shells use the real solid
%   harmonic polynomials of Ribaldone & Desmarais (2024, arXiv:2412.16733,
%   Table I -- a re-derivation of Schlegel & Frisch, Int. J. Quantum Chem.
%   54, 83 (1995)), combined from individually-normalized Cartesian
%   components with the overall normalization computed at runtime (not
%   hard-coded) from the exact combinatorial self-overlap of same-center
%   Cartesian Gaussians -- see PURE_HARMONIC_VALUE/CART_COMPONENT_VALUE
%   below), then combined via the MO coefficient column -- summed
%   directly into a single [Ngrid x 1] accumulator shell-by-shell, never
%   forming the full [Ngrid x Nbasis] basis matrix, to keep memory
%   bounded for large basis sets. The resulting scalar field is passed to
%   MATLAB's built-in ISOSURFACE/ISONORMALS/PATCH (no toolbox dependency).
%
%   S, P, SP, D, F, and G shells (pure or Cartesian) are supported; the
%   Cartesian shell component ordering (including the reversed order used
%   natively by Gaussian .fchk files for L>=4) follows IOData's
%   (iodata.formats.fchk) documented Gaussian-native convention. A basis
%   set containing H (or higher) shells raises a clear error rather than
%   silently omitting or mis-rendering those basis functions. F/G support
%   was validated on a real Au-containing basis (pure F and pure G on the
%   Au atom): self-normalization integrals of 0.998-1.000 for every MO
%   with significant weight on those shells.
%
%   A uniform rectangular grid under-resolves very tight (high-exponent)
%   primitives: deep core MOs, or high virtuals with significant weight
%   on tight core-region primitives, can show grid-dependent artifacts
%   near the nucleus even at fine GridSpacing (verified: self-
%   normalization integral non-monotonic/non-convergent for such MOs
%   under grid refinement, vs. clean 0.993-1.000 convergence for ordinary
%   valence/frontier MOs). This does not affect the chemically-relevant
%   use case of rendering valence/frontier (HOMO/LUMO-region) orbitals.
%
%   Automatic padding: some MOs (typically diffuse high virtuals) extend
%   well beyond the default 4 Bohr padding and would otherwise be visibly
%   clipped -- flat-cut isosurface faces at the edge of the plotted box.
%   After each isosurface extraction, the vertices are checked against
%   the grid boundary; if any lie within 1.5 grid cells of it, 'Padding'
%   is multiplied by 1.6 and the whole grid/evaluation/extraction is
%   redone, up to 3 times (subject to the same grid-size cap as a manual
%   'Padding'). If an explicit 'Padding' was given, this auto-expansion
%   is skipped and a warning recommends increasing it instead.
%
%   Example:
%       data = G16_fchk_read('4-NTP.fchk');
%       G_draw_mo_surface(data, 'HOMO');
%       G_draw_mo_surface(data, 'LUMO+1', 'IsoValue', 0.03, 'GridSpacing', 0.12);
%       G_draw_mo_surface(data, 'HOMO', 'Mode', 'contour');   % 2D nodal map
%       G_draw_mo_surface(data, 'HOMO', 'Mode', 'contour', 'PlaneOffset', 0.7);  % pi MO
%       G_draw_mo_surface(data, 'HOMO', 'SaveCube', 'homo.cube');   % also export
%       G_draw_mo_surface(data, 'HOMO', 'SurfaceStyle', 'grid', 'GridSpacing', 0.3);
%
%   See also G09_FCHK_READ, G16_FCHK_READ, G09_DRAW_MOLECULE, G16_DRAW_MOLECULE.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

% -------------------------------------------------------------------------
% Parse input arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'data');
addRequired(p,  'mo_index');
addParameter(p, 'Mode',         'surface',         @(s) any(strcmpi(s, {'surface','contour'})));
addParameter(p, 'Plane',        'auto',            @(s) any(strcmpi(s, {'auto','xy','xz','yz'})));
addParameter(p, 'PlaneOffset',  0,                 @isnumeric);
addParameter(p, 'IsoValue',     0.02,              @isnumeric);
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
parse(p, data, mo_index, varargin{:});

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
% 'FaceAlpha' always wins. 'grid' ignores FaceAlpha entirely (no face is
% drawn, only edges) -- see G_DRAW_DENSITY_SURFACE for the identical
% rationale.
if ~face_alpha_explicit && strcmp(surface_style, 'solid')
    face_alpha = 1.0;
end
if face_alpha_explicit && strcmp(surface_style, 'grid') && verbose
    fprintf('G_draw_mo_surface: ''FaceAlpha'' has no effect with ''SurfaceStyle'',''grid'' (no face is drawn, only mesh edges).\n');
end

if ~isempty(save_cube) && strcmp(mode, 'contour')
    error('G_draw_mo_surface: ''SaveCube'' needs the full 3D grid built in ''Mode'',''surface'' -- not available in ''Mode'',''contour'' (a single 2D plane).');
end

% -------------------------------------------------------------------------
% Validate input struct
% -------------------------------------------------------------------------
required = {'filename', 'mol', 'Nbasis', 'alpha_MO_coeff', 'xyz_bohr'};
for k = 1:numel(required)
    if ~isfield(data, required{k})
        error('G_draw_mo_surface: data must be the struct returned by G09_fchk_read/G16_fchk_read (missing "%s" field).', required{k});
    end
end
if ~isfile(data.filename)
    error('G_draw_mo_surface: data.filename (%s) not found -- this function re-reads the raw basis-set sections directly from the .fchk file.', data.filename);
end

aobasis = g_read_aobasis_from_fchk(data.filename);
if isempty(aobasis.shell_types)
    error('G_draw_mo_surface: no basis-set sections found in %s -- is this a valid .fchk file?', data.filename);
end

g_check_supported_shells(aobasis.shell_types);

homo_idx = [];
if isfield(data, 'HOMO_idx')
    homo_idx = data.HOMO_idx;
end
idx = resolve_mo_index(mo_index, homo_idx, data.Nbasis);

alpha_MO = reshape(data.alpha_MO_coeff, data.Nbasis, data.Nbasis);   % columns = MOs
mo_col   = alpha_MO(:, idx);

a0 = 0.529177210544;   % Bohr -> Angstrom (CODATA), consistent with the toolbox
xyz_bohr = data.xyz_bohr;

% -------------------------------------------------------------------------
% Set up figure/axes (shared by both modes)
% -------------------------------------------------------------------------
if isempty(ax)
    fig = figure('Color', 'white', 'Name', 'MO surface', 'NumberTitle', 'off');
    ax  = axes('Parent', fig);
end
hold(ax, 'on');

h.pos = [];
h.neg = [];
h.contour = [];
h.zero    = [];

if strcmp(mode, 'surface')
    % ---------------------------------------------------------------------
    % Build real-space grid (Bohr, for evaluation; Angstrom, for plotting)
    % and extract the isosurfaces, auto-expanding Padding and retrying if
    % the surface is clipped by the grid boundary -- unless the user gave
    % an explicit Padding, in which case their choice is respected and a
    % warning is issued instead. This avoids silently truncated diffuse
    % orbitals (e.g. high virtuals) without requiring manual Padding tuning.
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
                error(['G_draw_mo_surface: requested grid has %d points (%d x %d x %d), ' ...
                       'too large. Increase ''GridSpacing'' or decrease ''Padding''.'], ...
                      Ngrid, numel(gx), numel(gy), numel(gz));
            end
            warning('G_draw_mo_surface: stopped auto-expanding ''Padding'' at %.1f Bohr to stay under the %d-point grid limit; the surface may still be clipped. Increase ''GridSpacing'' to allow more room.', pad, max_grid);
            break
        end

        if verbose
            fprintf('G_draw_mo_surface: MO %d, grid %d x %d x %d = %d points (spacing=%.2f Bohr, padding=%.1f Bohr)\n', ...
                idx, numel(gx), numel(gy), numel(gz), Ngrid, spacing, pad);
        end

        [X, Y, Z] = meshgrid(gx, gy, gz);
        grid_pts_bohr = [X(:), Y(:), Z(:)];

        tic;
        mo_vals = g_eval_mo_on_grid(aobasis, mo_col, grid_pts_bohr);
        V = reshape(mo_vals, size(X));
        if verbose
            fprintf('  MO evaluation: %.1f s\n', toc);
        end

        Xa = X * a0;
        Ya = Y * a0;
        Za = Z * a0;

        fv_pos = isosurface(Xa, Ya, Za, V, isoval);
        fv_neg = isosurface(Xa, Ya, Za, V, -isoval);

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
            warning('G_draw_mo_surface: the isosurface appears to be clipped by the grid boundary (Padding=%.1f Bohr). Try a larger ''Padding''.', pad);
        else
            warning('G_draw_mo_surface: the isosurface still appears clipped after auto-expanding ''Padding'' to %.1f Bohr (%d attempts). This MO may be unusually diffuse; try an even larger ''Padding'' explicitly, or a larger ''GridSpacing'' to afford more room within the grid-size limit.', pad, attempt);
        end
    end

    if ~isempty(save_cube)
        title_line = sprintf('MO %d amplitude -- G_Utility export', idx);
        comment_line = sprintf('psi(r), 1/Bohr^1.5, grid spacing=%.3f Bohr', spacing);
        g_write_cube_file(save_cube, title_line, comment_line, data.mol.Z, xyz_bohr, gx, gy, gz, V, idx);
        if verbose
            fprintf('  Cube file written: %s\n', save_cube);
        end
    end

    % ---------------------------------------------------------------------
    % Draw isosurfaces
    % ---------------------------------------------------------------------
    if ~isempty(fv_pos.vertices)
        pos_style_args = surface_style_args(surface_style, pos_color, face_alpha);
        h.pos = patch(ax, fv_pos, pos_style_args{:}, ...
            'FaceLighting', 'gouraud', 'DisplayName', '+');
        isonormals(Xa, Ya, Za, V, h.pos);
    end

    if ~isempty(fv_neg.vertices)
        neg_style_args = surface_style_args(surface_style, neg_color, face_alpha);
        h.neg = patch(ax, fv_neg, neg_style_args{:}, ...
            'FaceLighting', 'gouraud', 'DisplayName', '-');
        isonormals(Xa, Ya, Za, V, h.neg);
    end

    if isempty(fv_pos.vertices) && isempty(fv_neg.vertices)
        warning('G_draw_mo_surface: no isosurface found at IsoValue=%.4g -- the MO amplitude may not reach this level anywhere on the grid. Try a smaller ''IsoValue'' or larger ''Padding''.', isoval);
    end

    % ---------------------------------------------------------------------
    % Overlay molecule (drawn after the surface so its own AXIS TIGHT call
    % fits both the surface and the atoms). Uses whichever of
    % G09_draw_molecule/G16_draw_molecule is found on the path -- both
    % accept the same data.mol struct.
    % ---------------------------------------------------------------------
    if show_mol
        if exist('G16_draw_molecule', 'file') == 2
            draw_mol_fcn = @G16_draw_molecule;
        elseif exist('G09_draw_molecule', 'file') == 2
            draw_mol_fcn = @G09_draw_molecule;
        else
            error('G_draw_mo_surface: ''ShowMolecule'' is true but neither G16_draw_molecule nor G09_draw_molecule is on the MATLAB path.');
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
    % Build the cutting plane, evaluate the MO on a 2D grid within it, and
    % render a filled amplitude contour with an explicit nodal (zero) line
    % -- see PLANE_BASIS/DRAW_2D_ATOMS below.
    % ---------------------------------------------------------------------
    [center, u, v, n] = g_plane_basis(xyz_bohr, plane_choice);
    offset_bohr = plane_offset / a0;

    atom_uv = (xyz_bohr - center) * [u, v];   % [Nat x 2], Bohr -- unaffected by
                                               % PlaneOffset: atoms are always
                                               % shown at their true in-plane
                                               % position for reference, even
                                               % when the MO is evaluated on a
                                               % parallel, offset plane.
    lo_uv = min(atom_uv, [], 1) - pad;
    hi_uv = max(atom_uv, [], 1) + pad;
    gs = lo_uv(1):spacing:hi_uv(1);
    gt = lo_uv(2):spacing:hi_uv(2);
    Ngrid = numel(gs) * numel(gt);

    max_grid = 4e6;
    if Ngrid > max_grid
        error(['G_draw_mo_surface: requested contour grid has %d points (%d x %d), ' ...
               'too large. Increase ''GridSpacing'' or decrease ''Padding''.'], ...
              Ngrid, numel(gs), numel(gt));
    end

    if verbose
        fprintf('G_draw_mo_surface: MO %d, contour grid %d x %d = %d points (spacing=%.2f Bohr, padding=%.1f Bohr, plane=%s, offset=%.2f %c)\n', ...
            idx, numel(gs), numel(gt), Ngrid, spacing, pad, plane_choice, plane_offset, 197);
    end

    [S, T] = meshgrid(gs, gt);
    grid_pts_bohr = (center + offset_bohr*n') + S(:)*u' + T(:)*v';

    tic;
    mo_vals = g_eval_mo_on_grid(aobasis, mo_col, grid_pts_bohr);
    V2 = reshape(mo_vals, size(S));
    if verbose
        fprintf('  MO evaluation: %.1f s\n', toc);
    end

    vmax = max(abs(V2(:)));

    if vmax < 1e-3
        warning(['G_draw_mo_surface: the MO amplitude is essentially zero (max %.2e a.u.) everywhere on ' ...
                 'this plane. A common cause: this is a pi-symmetry orbital and the plane (offset=%.2f %c) ' ...
                 'coincides with -- or sits very close to -- one of its nodal planes (e.g. the molecular ' ...
                 'plane itself, for a typical aromatic pi MO). Try a nonzero ''PlaneOffset'' (e.g. 0.5-1.0 ' ...
                 '%c) to probe a parallel plane where the pi density actually peaks.'], ...
                vmax, plane_offset, 197, 197);
    end

    Sa = S * a0;
    Ta = T * a0;

    % Diverging colormap (NegColor -> white -> PosColor), zero pinned to
    % white regardless of the actual data range (via symmetric CLim).
    n_cmap = 256;
    cmap = [interp1([0 1], [neg_color; 1 1 1], linspace(0,1,n_cmap/2)); ...
            interp1([0 1], [1 1 1; pos_color], linspace(0,1,n_cmap/2))];
    colormap(ax, cmap);

    if vmax == 0
        vmax = 1;   % degenerate (all-zero) field -- avoid a zero-width CLim
    end

    n_levels = 21;
    [~, h.contour] = contourf(ax, Sa, Ta, V2, linspace(-vmax, vmax, n_levels), 'LineStyle', 'none');
    caxis(ax, [-vmax, vmax]);   %#ok<CAXIS> -- kept over CLIM for MATLAB R2021b compatibility

    if vmax > 0 && any(V2(:) > 0, 'all') && any(V2(:) < 0, 'all')
        [~, h.zero] = contour(ax, Sa, Ta, V2, [0 0], 'Color', [0.15 0.15 0.15], 'LineWidth', 1.3);
    end

    cb = colorbar(ax);
    cb.Label.String = 'MO amplitude (a.u.)';

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
    tag = '';
    if ~isempty(homo_idx)
        if idx == homo_idx
            tag = ' (HOMO)';
        elseif idx == homo_idx + 1
            tag = ' (LUMO)';
        elseif idx < homo_idx
            tag = sprintf(' (HOMO-%d)', homo_idx - idx);
        else
            tag = sprintf(' (LUMO+%d)', idx - homo_idx - 1);
        end
    end
    e_str = '';
    if isfield(data, 'alpha_orb_energies') && idx <= numel(data.alpha_orb_energies)
        e_str = sprintf(', E = %.4f Ha', data.alpha_orb_energies(idx));
    end
    if strcmp(mode, 'surface')
        fig_title = sprintf('MO %d%s%s, iso=\\pm%.3g', idx, tag, e_str, isoval);
    else
        fig_title = sprintf('MO %d%s%s, plane=%s', idx, tag, e_str, plane_choice);
    end
end
title(ax, fig_title, 'Interpreter', 'tex', 'FontSize', 11);

hold(ax, 'off');

end % G_draw_mo_surface


% =========================================================================
%  Local functions
% =========================================================================

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

function idx = resolve_mo_index(mo_spec, homo_idx, n_basis)
%RESOLVE_MO_INDEX  Accepts a positive integer, or 'HOMO'/'LUMO'/'HOMO-n'/
%   'LUMO+n', and returns the 1-based MO column index.
    if isnumeric(mo_spec)
        idx = round(mo_spec);
    elseif ischar(mo_spec) || isstring(mo_spec)
        s = upper(strtrim(char(mo_spec)));
        tok = regexp(s, '^(HOMO|LUMO)([+-]\d+)?$', 'tokens', 'once');
        if isempty(tok)
            error('G_draw_mo_surface: mo_index string must be numeric, or ''HOMO''/''LUMO''/''HOMO-n''/''LUMO+n''.');
        end
        if isempty(homo_idx)
            error('G_draw_mo_surface: mo_index=''%s'' requires data.HOMO_idx, which is missing.', s);
        end
        base = tok{1};
        offs = tok{2};
        if isempty(offs)
            n = 0;
        else
            n = str2double(offs);
        end
        if strcmp(base, 'HOMO')
            idx = homo_idx + n;
        else
            idx = homo_idx + 1 + n;
        end
    else
        error('G_draw_mo_surface: mo_index must be a positive integer or a string like ''HOMO'', ''LUMO'', ''HOMO-1'', ''LUMO+2''.');
    end
    if idx < 1 || idx > n_basis
        error('G_draw_mo_surface: resolved MO index %d is out of range [1, %d].', idx, n_basis);
    end
end

