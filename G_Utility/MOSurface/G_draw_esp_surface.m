function h = G_draw_esp_surface(data, varargin)
% G_DRAW_ESP_SURFACE  Renders the molecular electrostatic potential (ESP)
%                      mapped onto the total-electron-density isosurface
%                      -- a G_Utility equivalent of GaussView's classic
%                      "ESP mapped on electron density" surface.
%
%   G_draw_esp_surface(data)
%   G_draw_esp_surface(data, 'Name', Value, ...)
%   G_draw_esp_surface(data, 'CompareTo', data2, ...)
%   h = G_draw_esp_surface(...)
%
%   This is a standalone G_Utility companion function (see
%   G_DRAW_MO_SURFACE for the shared design rationale) that reuses
%   G_DRAW_DENSITY_SURFACE's machinery to build the density envelope,
%   then colours it by the electrostatic potential
%       V(r) = sum_A Z_A/|r-R_A|  -  integral[ rho(r')/|r-r'| dr' ]
%   (nuclear minus electronic contribution; standard sign convention --
%   V<0 = electron-rich/nucleophilic, V>0 = electron-poor/electrophilic).
%   Z_A is read from the .fchk "Nuclear charges" section (NOT the atomic
%   number): for an ECP (effective core potential) atom, this is the
%   reduced effective charge seen by the explicit electrons, not the full
%   atomic number -- using the atomic number instead would silently give
%   a badly wrong nuclear term for any ECP-treated heavy atom (verified
%   on a real Au-containing test file: "Nuclear charges" = 19 for Au,
%   vs. atomic number 79 -- a 60-electron core ECP).
%
%   The electronic term is evaluated by the EXACT analytic
%   McMurchie-Davidson Gaussian integral scheme (G_EVAL_ESP_ANALYTIC),
%   validated end-to-end against a real Gaussian-generated ESP cube file
%   (correlation 1.0000000, RMS error ~0% of the true signal span at the
%   actual density-surface distance range) -- see the "Method" and
%   "Accuracy" notes below. G_EVAL_ESP_ANALYTIC only supports S, P, SP,
%   and pure D shells; for molecules with Cartesian D or F/G shells (or
%   if 'ESPMethod' is explicitly set to 'numeric'), this function falls
%   back to the older NUMERICAL Coulomb-grid-sum method, which remains
%   available but was found to have substantial error (~54% of the true
%   ESP signal span) at the realistic density-surface distance -- treat
%   its output as qualitative only.
%
%   Input:
%       data      - struct returned by G09_FCHK_READ or G16_FCHK_READ.
%                   Required fields: .filename, .mol, .Nbasis,
%                   .alpha_MO_coeff, .xyz_bohr, .Nalpha, .Nbeta, .Nelec.
%                   Closed-shell only (.Nalpha == .Nbeta), as in
%                   G_DRAW_DENSITY_SURFACE.
%
%   Optional parameters (Name-Value):
%       'CompareTo'     - a second data struct (same requirements as
%                         DATA). If given, colours the surface by the ESP
%                         DIFFERENCE V(CompareTo) - V(data) instead of the
%                         absolute ESP of data alone -- e.g. field-on
%                         minus field-off, or two different geometries/
%                         electronic states of the same molecule, to see
%                         where the electrostatic potential increases/
%                         decreases. Both ESPs are evaluated independently
%                         (each with its own analytic/numerical fallback,
%                         nuclear charges, and density matrix) at the SAME
%                         vertices -- built from DATA's density surface
%                         shape only, exactly as G_DRAW_DENSITY_SURFACE's
%                         'CompareTo' -- so the comparison is only
%                         physically meaningful if data and CompareTo
%                         share the same atomic coordinate frame (default:
%                         [], absolute ESP of data only)
%       'IsoValue'      - density isovalue defining the surface shape, in
%                         electrons/Bohr^3 (default: 0.001, the classic
%                         Bader molecular envelope, as in
%                         G_DRAW_DENSITY_SURFACE)
%       'GridSpacing'   - real-space grid spacing for the surface SHAPE,
%                         in Bohr (default: 0.15)
%       'Padding'       - grid padding for the surface SHAPE, in Bohr
%                         (default: 4.0; auto-expanded if clipped, as in
%                         G_DRAW_MO_SURFACE)
%       'MaxVertices'   - the density isosurface mesh is decimated (via
%                         MATLAB's built-in REDUCEPATCH) to at most this
%                         many vertices before computing the ESP at each
%                         one, since ESP is smooth/slowly-varying and a
%                         full-resolution mesh (often >10000 vertices) is
%                         unnecessary and would make the O(vertices x
%                         ESP-grid-points) Coulomb sum far slower than
%                         needed for a visually smooth colour map
%                         (default: 3000)
%       'ESPMethod'     - 'auto' (default): use the exact analytic
%                         McMurchie-Davidson method (G_EVAL_ESP_ANALYTIC)
%                         whenever the basis is within its supported
%                         scope (S/P/SP/pure-D), otherwise fall back to
%                         the numerical Coulomb-grid-sum method with a
%                         notice; 'analytic': force the analytic method,
%                         erroring out (via G_EVAL_ESP_ANALYTIC) if the
%                         basis has unsupported shells; 'numeric': force
%                         the older, faster but substantially less
%                         accurate Coulomb-grid-sum method (see
%                         "Accuracy" below) regardless of basis.
%       'ESPGridSpacing'- ONLY used by the numerical (grid-sum) method,
%                         i.e. when 'ESPMethod' is 'numeric' or an 'auto'
%                         fallback was triggered: spacing of the SEPARATE
%                         density grid used for the electronic Coulomb
%                         sum, in Bohr (default: 0.20). An earlier design
%                         assumption that this could safely be much
%                         coarser than 'GridSpacing', since Coulomb
%                         potential is a smooth, long-range-integrated
%                         quantity, turned out empirically WRONG -- see
%                         "Accuracy" below; 0.20 is the finest value
%                         tested (this function does not itself get any
%                         more expensive for finer spacing beyond the
%                         Coulomb-grid density evaluation cost, so a
%                         still finer value is a reasonable thing to try
%                         if runtime allows)
%       'ESPPadding'    - ONLY used by the numerical (grid-sum) method
%                         (see 'ESPGridSpacing' above): padding for that
%                         same grid, in Bohr (default: 6.0). Empirically,
%                         extending this further makes little difference
%                         -- verified on a real molecule, doubling
%                         'ESPPadding' left the grid's captured electron
%                         count essentially unchanged, since the density
%                         already falls to near-zero well within the
%                         default padding; 'ESPGridSpacing' is what
%                         matters (see "Accuracy" below)
%       'DistanceFloor' - ONLY used by the numerical (grid-sum) method
%                         (see 'ESPGridSpacing' above): minimum distance
%                         (Bohr) used when dividing by |r-r'| in the
%                         Coulomb sum, to avoid a spurious near-singular
%                         spike when an ESP-grid density point happens to
%                         land very close to a surface vertex (default:
%                         [] = 0.5*'ESPGridSpacing')
%       'SaveCube'      - filename; if given, ALSO evaluates ESP on a
%                         full, independent volumetric grid (see
%                         'CubeSpacing'/'CubePadding' below, NOT the
%                         surface-shape grid, and NOT just the decimated
%                         surface vertices) and writes it to a
%                         Gaussian-format .cube file. This is a genuinely
%                         separate, additional computation -- unlike
%                         G_DRAW_DENSITY_SURFACE/G_DRAW_MO_SURFACE, ESP is
%                         only ever evaluated at the (few thousand at
%                         most) decimated surface vertices, never on a
%                         full grid, so there is no existing volumetric
%                         array to reuse for free. Runtime scales with the
%                         cube grid's point count and (for the default
%                         analytic method) roughly with Nbasis^2 per
%                         point -- can be genuinely slow for a large basis
%                         (default: '', no file written)
%       'CubeSpacing'   - grid spacing for 'SaveCube', in Bohr. Coarser
%                         than 'GridSpacing'/'ESPGridSpacing' by default,
%                         since a fine cube (matching, say, 0.15-0.20 Bohr)
%                         can mean hundreds of thousands of points, each
%                         costing a full shell-pair loop with the
%                         analytic method (default: 0.30)
%       'CubePadding'   - padding for 'SaveCube', in Bohr (default: 4.0)
%       'PosColor'      - colour for positive (electron-poor,
%                         electrophilic) ESP (default: [0.10 0.40 0.85],
%                         blue -- the standard ESP convention)
%       'NegColor'      - colour for negative (electron-rich,
%                         nucleophilic) ESP (default: [0.85 0.15 0.10],
%                         red -- the standard ESP convention)
%       'FaceAlpha'     - surface opacity, 0-1. Its own default (1.0) is
%                         used as-is with 'SurfaceStyle','solid' (the
%                         default here, unlike G_DRAW_MO_SURFACE/
%                         G_DRAW_DENSITY_SURFACE's 'transparent' default
%                         -- an ESP map is conventionally an opaque
%                         painted surface); overridden to 0.55 with
%                         'SurfaceStyle','transparent' UNLESS given
%                         explicitly (which always wins); has no effect
%                         with 'SurfaceStyle','grid'. EXCEPT: if
%                         'ShowMolecule' is true, 'SurfaceStyle' is
%                         'solid' (the default), and 'FaceAlpha' is left
%                         at its default, it is automatically lowered to
%                         0.6 instead, since at FaceAlpha=1 the CPK model
%                         would be completely invisible -- entirely
%                         enclosed by the opaque surface -- defeating the
%                         purpose of asking for it. Pass 'FaceAlpha'
%                         explicitly to keep it opaque regardless. (No
%                         such override is needed for 'transparent'
%                         -- already translucent -- or 'grid' -- no face
%                         to hide the molecule behind in the first
%                         place.)
%       'SurfaceStyle'  - 'solid' (default) | 'transparent' | 'grid' --
%                         as in GaussView's own surface-style choice.
%                         'grid' draws the surface as a wireframe (mesh
%                         edges only, no face fill), coloured by
%                         interpolating the same ESP colormap along the
%                         edges (EdgeColor,'interp') rather than a flat
%                         colour, since ESP (unlike the density/MO
%                         lobes) has no single per-region colour to
%                         begin with. Neither the surface-shape mesh nor
%                         the decimated ('MaxVertices') vertex set
%                         changes with 'SurfaceStyle' -- 'grid' can
%                         still be dense; reduce 'MaxVertices' for a
%                         more readable wireframe
%       'ShowMolecule'  - overlay the CPK ball-and-stick model beneath
%                         the ESP surface (default: false). See
%                         'FaceAlpha' above: requesting this without also
%                         specifying an explicit 'FaceAlpha' switches a
%                         'SurfaceStyle','solid' surface to translucent
%                         so the molecule is actually visible.
%       'ShowLabels'    - as in G_DRAW_MO_SURFACE (default: false)
%       'AtomScale'     - as in G_DRAW_MO_SURFACE (default: 0.35)
%       'Title'         - figure title (default: auto)
%       'Ax'            - existing axes handle (default: new figure)
%       'Verbose'       - print grid size and timing info, including
%                         progress through the Coulomb-sum batches, since
%                         this is by far the slowest function in this
%                         toolbox -- roughly half a minute for a ~15-atom
%                         molecule at the default settings (default: true)
%
%   Output:
%       h - struct with field .surf (the coloured surface patch handle).
%
%   Method: (1) build the density isosurface exactly as in
%   G_DRAW_DENSITY_SURFACE, decimate its mesh to 'MaxVertices'; (2) at
%   each (decimated) surface vertex, V_nuclear is the exact analytic
%   point-charge sum over nuclei (Z_A/|r-R_A|, from the .fchk "Nuclear
%   charges"), and V_electronic is evaluated by G_EVAL_ESP_ANALYTIC
%   (exact McMurchie-Davidson Gaussian integrals against the density
%   matrix) when the basis is within its supported scope, or else by the
%   numerical Coulomb-grid-sum fallback described under "Accuracy"
%   below; (3) V = V_nuclear - V_electronic is painted onto the mesh via
%   per-vertex FaceVertexCData.
%
%   Accuracy: the default analytic method (G_EVAL_ESP_ANALYTIC) is exact
%   up to double-precision/basis-set truncation -- validated end-to-end
%   against a real Gaussian-generated ESP cube file at the actual
%   density-isosurface distance range (2.4-4.4 Bohr from the nearest
%   nucleus): correlation 1.0000000, RMS error ~8.7e-7 Hartree/e (0.0% of
%   the true ~0.14 Hartree/e signal span). It only supports S, P, SP, and
%   pure D shells, however; for a basis with Cartesian D or F/G shells
%   this function automatically falls back (with a printed notice, if
%   'Verbose') to the older NUMERICAL Coulomb-grid-sum method described
%   below, or errors out if 'ESPMethod' is explicitly 'analytic'.
%
%   The numerical fallback evaluates the density once on a separate,
%   coarser, wider grid via G_EVAL_DENSITY_ON_GRID, then UNIFORMLY
%   RESCALES it so its total integral exactly equals the known-exact
%   electron count data.Nelec (not optional bookkeeping -- it corrects a
%   real, sometimes large systematic bias), then sums
%   rho(r'_k) dV / max(|r-r'_k|, DistanceFloor) over the (rescaled) grid
%   in vertex batches to bound memory. This numerical method was found,
%   by the same real-cube-file validation described above, to have
%   substantial error at the realistic surface distance (correlation
%   0.03, RMS error 54% of the true signal span) -- adequate at most for
%   the qualitative electrophilic/nucleophilic pattern an ESP map is
%   normally read for, NOT for quantitative values. A uniform grid
%   systematically MIS-ESTIMATES the raw density right at each nuclear
%   cusp: verified on a real molecule, the total electron count captured
%   by a naive grid sum ranged from -2% to +31% depending on
%   'ESPGridSpacing' alone, NON-monotonically (refining the grid does not
%   straightforwardly reduce this -- the classic symptom of aliasing a
%   sharply peaked integrand on a uniform grid, the same phenomenon
%   documented for individual MOs in G_DRAW_MO_SURFACE); the rescaling
%   step removes the dominant, roughly-uniform component of this error,
%   but not its remaining spatial (atom-to-atom) variation. The remaining
%   spacing sensitivity is not just cosmetic: on the same real molecule,
%   'ESPGridSpacing' 0.5 (fast) gave a rescaled ESP ranging +-1.5 to +-3
%   Hartree/e -- unphysically large for a vdW-like density envelope, real
%   organic molecules are typically +-0.05 to +-0.15 -- and only showed a
%   broad, smeared +/- split; 'ESPGridSpacing' 0.20 (the default, ~10x
%   slower) gave the physically-reasonable +-0.14 range and clearly
%   resolved localized features expected from the underlying chemistry
%   (e.g. a distinct negative well right at a nitro group's
%   lone-pair-bearing oxygens, verified on the same test molecule) that
%   the coarser setting smeared away entirely. 'DistanceFloor' separately
%   guards against a spurious spike when a Coulomb-grid point happens to
%   land very close to a surface vertex.
%
%   S, P, SP, D, F, and G shells (pure or Cartesian) are supported for
%   the density surface SHAPE itself, as in G_DRAW_DENSITY_SURFACE --
%   only the analytic ESP method (see above) has the narrower S/P/SP/
%   pure-D scope.
%
%   Example:
%       data = G16_fchk_read('4-NTP.fchk');
%       G_draw_esp_surface(data);
%       G_draw_esp_surface(data, 'MaxVertices', 1500);  % faster (fewer Coulomb-sum
%                                                        % vertices), same grid accuracy
%       data2 = G16_fchk_read('4-NTP_reduced.fchk');
%       G_draw_esp_surface(data, 'CompareTo', data2);   % ESP difference map
%       G_draw_esp_surface(data, 'SaveCube', 'esp.cube');           % also export
%       G_draw_esp_surface(data, 'SaveCube', 'esp.cube', 'CubeSpacing', 0.5);  % faster/coarser
%       G_draw_esp_surface(data, 'SurfaceStyle', 'grid', 'MaxVertices', 800);
%
%   See also G_DRAW_DENSITY_SURFACE, G_DRAW_MO_SURFACE, G09_FCHK_READ, G16_FCHK_READ.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

% -------------------------------------------------------------------------
% Parse input arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'data');
addParameter(p, 'IsoValue',      0.001,             @isnumeric);
addParameter(p, 'GridSpacing',   0.15,              @isnumeric);
addParameter(p, 'Padding',       4.0,               @isnumeric);
addParameter(p, 'MaxVertices',   3000,              @isnumeric);
addParameter(p, 'CompareTo',     [],                @(s) isempty(s) || isstruct(s));
addParameter(p, 'ESPMethod',     'auto',            @(x) ischar(x) && ismember(lower(x), {'auto','analytic','numeric'}));
addParameter(p, 'ESPGridSpacing',0.20,              @isnumeric);
addParameter(p, 'ESPPadding',    6.0,               @isnumeric);
addParameter(p, 'DistanceFloor', [],                @(x) isempty(x) || isnumeric(x));
addParameter(p, 'SaveCube',      '',                @ischar);
addParameter(p, 'CubeSpacing',   0.30,              @isnumeric);
addParameter(p, 'CubePadding',   4.0,               @isnumeric);
addParameter(p, 'PosColor',      [0.10 0.40 0.85],  @isnumeric);
addParameter(p, 'NegColor',      [0.85 0.15 0.10],  @isnumeric);
addParameter(p, 'FaceAlpha',     1.0,               @isnumeric);
addParameter(p, 'SurfaceStyle',  'solid',           @(s) any(strcmpi(s, {'solid','transparent','grid'})));
addParameter(p, 'ShowMolecule',  false,             @islogical);
addParameter(p, 'ShowLabels',    false,             @islogical);
addParameter(p, 'AtomScale',     0.35,              @isnumeric);
addParameter(p, 'Title',         '',                @ischar);
addParameter(p, 'Ax',            [],                @ishandle);
addParameter(p, 'Verbose',       true,              @islogical);
parse(p, data, varargin{:});

isoval        = abs(p.Results.IsoValue);
spacing       = p.Results.GridSpacing;
pad           = p.Results.Padding;
max_vertices  = p.Results.MaxVertices;
data2         = p.Results.CompareTo;
compare_mode  = ~isempty(data2);
esp_method    = lower(p.Results.ESPMethod);
esp_spacing   = p.Results.ESPGridSpacing;
esp_pad       = p.Results.ESPPadding;
dist_floor    = p.Results.DistanceFloor;
save_cube     = p.Results.SaveCube;
cube_spacing  = p.Results.CubeSpacing;
cube_pad      = p.Results.CubePadding;
pos_color     = p.Results.PosColor;
neg_color     = p.Results.NegColor;
face_alpha    = p.Results.FaceAlpha;
surface_style = lower(p.Results.SurfaceStyle);
show_mol      = p.Results.ShowMolecule;
show_labels   = p.Results.ShowLabels;
atom_scale    = p.Results.AtomScale;
fig_title     = p.Results.Title;
ax            = p.Results.Ax;
verbose       = p.Results.Verbose;
padding_explicit    = ~ismember('Padding', p.UsingDefaults);
face_alpha_explicit = ~ismember('FaceAlpha', p.UsingDefaults);

% 'SurfaceStyle' only sets FaceAlpha's DEFAULT: 'solid' keeps the
% existing opaque default (1.0, unchanged); 'transparent' uses a
% translucent default (0.55, matching G_DRAW_DENSITY_SURFACE/
% G_DRAW_MO_SURFACE) unless FaceAlpha was given explicitly (which always
% wins); 'grid' ignores FaceAlpha entirely (no face is drawn, only mesh
% edges, coloured -- like the face would have been -- by interpolated
% ESP value; see 'grid' handling below).
if ~face_alpha_explicit && strcmp(surface_style, 'transparent')
    face_alpha = 0.55;
end
if face_alpha_explicit && strcmp(surface_style, 'grid') && verbose
    fprintf('G_draw_esp_surface: ''FaceAlpha'' has no effect with ''SurfaceStyle'',''grid'' (no face is drawn, only mesh edges).\n');
end

if isempty(dist_floor)
    dist_floor = 0.5 * esp_spacing;
end

% At the default FaceAlpha=1 (opaque, 'SurfaceStyle','solid'), a molecule
% overlay requested via 'ShowMolecule' is entirely invisible -- the CPK
% spheres/bonds are drawn (verified: G16_draw_molecule's legend does
% appear), just completely hidden inside/behind the opaque ESP surface
% that encloses them. Since an explicit 'ShowMolecule',true is a direct
% request to see the molecule, silently doing nothing about that
% contradiction is worse than picking a sensible translucent default --
% unless the user also gave an explicit 'FaceAlpha' (respected as-is), or
% chose 'SurfaceStyle','grid' (no face to hide the molecule behind in the
% first place, so this override does not apply there).
if show_mol && ~strcmp(surface_style, 'grid')
    if ~face_alpha_explicit && strcmp(surface_style, 'solid')
        face_alpha = 0.6;
        if verbose
            fprintf('G_draw_esp_surface: ''ShowMolecule'' is true, so ''FaceAlpha'' defaults to %.2g (not the usual opaque 1.0) so the molecule underneath is actually visible. Pass ''FaceAlpha'' explicitly to override.\n', face_alpha);
        end
    elseif face_alpha_explicit && face_alpha > 0.9
        warning('G_draw_esp_surface: ''ShowMolecule'' is true but ''FaceAlpha''=%.2g is nearly opaque -- the molecule will be mostly or entirely hidden underneath the surface.', face_alpha);
    end
end

% -------------------------------------------------------------------------
% Validate input struct(s) and prepare each one's occupied-MO coefficient
% matrix, density matrix, and nuclear charges (closed-shell only)
% -------------------------------------------------------------------------
[aobasis, occ_coeff, P_density, nuclear_charges] = validate_and_prepare(data, 'data');
if compare_mode
    [aobasis2, occ_coeff2, P_density2, nuclear_charges2] = validate_and_prepare(data2, 'CompareTo');
end

a0 = 0.529177210544;   % Bohr -> Angstrom (CODATA), consistent with the toolbox
xyz_bohr = data.xyz_bohr;

% -------------------------------------------------------------------------
% Build the density isosurface (shape only) -- same auto-expanding grid
% logic as G_draw_density_surface's 'surface' mode, single (positive)
% lobe only since density is never negative.
% -------------------------------------------------------------------------
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
            error(['G_draw_esp_surface: requested surface-shape grid has %d points (%d x %d x %d), ' ...
                   'too large. Increase ''GridSpacing'' or decrease ''Padding''.'], ...
                  Ngrid, numel(gx), numel(gy), numel(gz));
        end
        warning('G_draw_esp_surface: stopped auto-expanding ''Padding'' at %.1f Bohr to stay under the %d-point grid limit.', pad, max_grid);
        break
    end

    if verbose
        fprintf('G_draw_esp_surface: surface-shape grid %d x %d x %d = %d points (spacing=%.2f Bohr, padding=%.1f Bohr)\n', ...
            numel(gx), numel(gy), numel(gz), Ngrid, spacing, pad);
    end

    [X, Y, Z] = meshgrid(gx, gy, gz);
    grid_pts_bohr = [X(:), Y(:), Z(:)];

    tic;
    V = reshape(g_eval_density_on_grid(aobasis, occ_coeff, grid_pts_bohr), size(X));
    if verbose
        fprintf('  Density evaluation: %.1f s\n', toc);
    end

    Xa = X * a0; Ya = Y * a0; Za = Z * a0;
    fv = isosurface(Xa, Ya, Za, V, isoval);

    clipped = g_is_clipped(fv, Xa, Ya, Za, 1.5*spacing*a0);

    if ~clipped || padding_explicit || attempt >= max_attempts
        break
    end
    if verbose
        fprintf('  Isosurface touches the grid boundary -- expanding Padding %.1f -> %.1f Bohr and retrying.\n', pad, pad*pad_growth);
    end
    pad = pad * pad_growth;
end

if clipped
    warning('G_draw_esp_surface: the density isosurface appears clipped by the grid boundary (Padding=%.1f Bohr). Try a larger ''Padding''.', pad);
end

if isempty(fv.vertices)
    error('G_draw_esp_surface: no density isosurface found at IsoValue=%.4g -- try a smaller ''IsoValue'' or larger ''Padding''.', isoval);
end

% -------------------------------------------------------------------------
% Decimate the mesh (ESP is smooth; a full-resolution density mesh has
% far more vertices than needed for a visually smooth colour map, and
% every extra vertex costs one more Coulomb-sum evaluation).
% -------------------------------------------------------------------------
n_verts_full = size(fv.vertices, 1);
if n_verts_full > max_vertices
    reduce_frac = max_vertices / n_verts_full;
    fv = reducepatch(fv, reduce_frac);
    if verbose
        fprintf('  Mesh decimated: %d -> %d vertices\n', n_verts_full, size(fv.vertices,1));
    end
end
verts_ang  = fv.vertices;          % [Nv x 3], Angstrom
verts_bohr = verts_ang / a0;       % [Nv x 3], Bohr

% -------------------------------------------------------------------------
% ESP at each surface vertex -- exact analytic McMurchie-Davidson
% evaluation (G_EVAL_ESP_ANALYTIC) by default, falling back to the
% numerical Coulomb-grid-sum method (see "Accuracy" in the help text
% above) for basis sets outside the analytic method's S/P/SP/pure-D
% scope, or if 'ESPMethod' explicitly requests 'numeric'. In 'CompareTo'
% mode, both calculations are evaluated independently at the SAME
% vertices and subtracted.
% -------------------------------------------------------------------------
V_esp1 = compute_esp_at_points(aobasis, occ_coeff, P_density, nuclear_charges, ...
    xyz_bohr, data.mol.Natoms, data.Nelec, verts_bohr, esp_method, esp_spacing, esp_pad, dist_floor, verbose, 'data');

if compare_mode
    V_esp2 = compute_esp_at_points(aobasis2, occ_coeff2, P_density2, nuclear_charges2, ...
        data2.xyz_bohr, data2.mol.Natoms, data2.Nelec, verts_bohr, esp_method, esp_spacing, esp_pad, dist_floor, verbose, 'CompareTo');
    V_esp = V_esp2 - V_esp1;   % Hartree/e, CompareTo minus data
else
    V_esp = V_esp1;            % Hartree/e, standard sign convention
end

% -------------------------------------------------------------------------
% Optional: ESP on a full volumetric grid, written to a .cube file. Unlike
% the density/MO cube export (which reuses an already-computed grid for
% free), ESP is normally only ever evaluated at the decimated surface
% vertices above -- this is a genuinely separate, additional evaluation
% over a dedicated (by default coarser) grid.
% -------------------------------------------------------------------------
if ~isempty(save_cube)
    max_cube_grid = 2e6;   % much smaller than the 30e6 shape-grid cap:
                            % ESP is orders of magnitude more expensive
                            % per point than density/MO point evaluation
    lo_c = min(xyz_bohr, [], 1) - cube_pad;
    hi_c = max(xyz_bohr, [], 1) + cube_pad;
    gxc = lo_c(1):cube_spacing:hi_c(1);
    gyc = lo_c(2):cube_spacing:hi_c(2);
    gzc = lo_c(3):cube_spacing:hi_c(3);
    Ngrid_cube = numel(gxc) * numel(gyc) * numel(gzc);
    if Ngrid_cube > max_cube_grid
        error(['G_draw_esp_surface: ''SaveCube'' grid has %d points (%d x %d x %d), ' ...
               'too large for a per-point ESP evaluation (cap: %d). Increase ' ...
               '''CubeSpacing'' or decrease ''CubePadding''.'], ...
              Ngrid_cube, numel(gxc), numel(gyc), numel(gzc), max_cube_grid);
    end
    if verbose
        fprintf('G_draw_esp_surface: ''SaveCube'' grid %d x %d x %d = %d points (spacing=%.2f Bohr, padding=%.1f Bohr) -- this is a separate, additional evaluation and may take a while.\n', ...
            numel(gxc), numel(gyc), numel(gzc), Ngrid_cube, cube_spacing, cube_pad);
    end

    [Xc, Yc, Zc] = meshgrid(gxc, gyc, gzc);
    cube_pts_bohr = [Xc(:), Yc(:), Zc(:)];

    Vc1 = compute_esp_at_points(aobasis, occ_coeff, P_density, nuclear_charges, ...
        xyz_bohr, data.mol.Natoms, data.Nelec, cube_pts_bohr, esp_method, esp_spacing, esp_pad, dist_floor, verbose, 'data, cube grid');
    if compare_mode
        Vc2 = compute_esp_at_points(aobasis2, occ_coeff2, P_density2, nuclear_charges2, ...
            data2.xyz_bohr, data2.mol.Natoms, data2.Nelec, cube_pts_bohr, esp_method, esp_spacing, esp_pad, dist_floor, verbose, 'CompareTo, cube grid');
        Vcube = reshape(Vc2 - Vc1, size(Xc));
        title_line = 'ESP difference (CompareTo - data) -- G_Utility export';
    else
        Vcube = reshape(Vc1, size(Xc));
        title_line = 'Electrostatic potential -- G_Utility export';
    end
    comment_line = sprintf('V(r), Hartree/e, grid spacing=%.3f Bohr', cube_spacing);
    g_write_cube_file(save_cube, title_line, comment_line, data.mol.Z, xyz_bohr, gxc, gyc, gzc, Vcube);
    if verbose
        fprintf('  Cube file written: %s\n', save_cube);
    end
end

% -------------------------------------------------------------------------
% Set up figure/axes and draw the coloured surface
% -------------------------------------------------------------------------
if isempty(ax)
    fig = figure('Color', 'white', 'Name', 'ESP surface', 'NumberTitle', 'off');
    ax  = axes('Parent', fig);
end
hold(ax, 'on');

n_cmap = 256;
cmap = [interp1([0 1], [neg_color; 1 1 1], linspace(0,1,n_cmap/2)); ...
        interp1([0 1], [1 1 1; pos_color], linspace(0,1,n_cmap/2))];
colormap(ax, cmap);
vmax = g_pctile_local(abs(V_esp), 99);
if vmax == 0, vmax = 1; end

% 'grid' colours the wireframe edges by the same interpolated ESP value
% the faces would have used (EdgeColor,'interp' with FaceVertexCData is
% supported directly by PATCH, exactly like FaceColor,'interp') --
% unlike G_DRAW_DENSITY_SURFACE/G_DRAW_MO_SURFACE's flat PosColor/
% NegColor wireframe, since ESP has no single per-lobe colour to begin
% with.
if strcmp(surface_style, 'grid')
    esp_style_args = {'FaceColor', 'none', 'EdgeColor', 'interp', 'FaceAlpha', 1};
else   % 'solid' / 'transparent'
    esp_style_args = {'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', face_alpha};
end

h.surf = patch(ax, 'Vertices', verts_ang, 'Faces', fv.faces, ...
    'FaceVertexCData', V_esp, esp_style_args{:}, 'FaceLighting', 'gouraud', ...
    'VertexNormalsMode', 'auto');   % MATLAB computes shading normals from
                                     % face geometry; no source volume is
                                     % available post-REDUCEPATCH for ISONORMALS
caxis(ax, [-vmax, vmax]);   %#ok<CAXIS> -- kept over CLIM for MATLAB R2021b compatibility

cb = colorbar(ax);
if compare_mode
    cb.Label.String = 'ESP difference, CompareTo - data (Hartree/e)';
else
    cb.Label.String = 'Electrostatic potential (Hartree/e)';
end

% Lighting/view setup delegated to G09/G16_draw_molecule when
% 'ShowMolecule' is true, to avoid setting up two independent pairs of
% CAMLIGHTs (this function's own, plus draw_molecule's) -- which would
% over-brighten the surface relative to G_draw_mo_surface/
% G_draw_density_surface (see their own identical if/else pattern).
if show_mol
    if exist('G16_draw_molecule', 'file') == 2
        draw_mol_fcn = @G16_draw_molecule;
    elseif exist('G09_draw_molecule', 'file') == 2
        draw_mol_fcn = @G09_draw_molecule;
    else
        error('G_draw_esp_surface: ''ShowMolecule'' is true but neither G16_draw_molecule nor G09_draw_molecule is on the MATLAB path.');
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

if isempty(fig_title)
    if compare_mode
        fig_title = sprintf('ESP difference on density surface (iso=%.3g)', isoval);
    else
        fig_title = sprintf('ESP on density surface (iso=%.3g)', isoval);
    end
end
title(ax, fig_title, 'Interpreter', 'tex', 'FontSize', 11);

hold(ax, 'off');

end % G_draw_esp_surface


% =========================================================================
%  Local functions
% =========================================================================

function Z = read_nuclear_charges(filename, Natoms)
%READ_NUCLEAR_CHARGES  Reads the .fchk "Nuclear charges" section directly
%   (a minimal, single-purpose parser, in the same self-contained spirit
%   as G_READ_AOBASIS_FROM_FCHK) -- NOT the same as "Atomic numbers": for
%   an ECP (effective core potential) atom, "Nuclear charges" gives the
%   reduced effective charge seen by the explicit electrons, which is the
%   value the electrostatic potential's nuclear term must use to be
%   correct (verified on a real Au-ECP test file: Nuclear charges = 19
%   for Au, vs. atomic number 79).
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
        if ~strcmpi(strtrim(tok{1}), 'Nuclear charges'), continue; end
        nvals = str2double(tok{3});
        vals = [];
        k2 = k + 1;
        while numel(vals) < nvals && k2 <= N
            ln2 = strtrim(lines{k2});
            if isempty(ln2), k2 = k2+1; continue; end
            if numel(ln2) > 44 && ~isempty(regexp(ln2, header_re, 'once'))
                break
            end
            vals = [vals; sscanf(ln2, '%f')]; %#ok<AGROW>
            k2 = k2 + 1;
        end
        Z = vals(1:min(end,nvals));
        Z = Z(:);
        if numel(Z) ~= Natoms
            error('G_draw_esp_surface: "Nuclear charges" section in %s has %d entries, expected %d (Natoms).', filename, numel(Z), Natoms);
        end
        return
    end
    error('G_draw_esp_surface: "Nuclear charges" section not found in %s -- is this a valid .fchk file?', filename);
end

function [aobasis, occ_coeff, P_density, nuclear_charges] = validate_and_prepare(data, argname)
%VALIDATE_AND_PREPARE  Validates one data struct (either the primary DATA
%   or a 'CompareTo' struct), re-reads its raw basis-set data, checks it
%   is closed-shell, and returns everything G_EVAL_ESP_ANALYTIC/the
%   numerical fallback need: the AO basis, the occupied-MO coefficient
%   matrix [Nbasis x Nalpha], the density matrix P = 2*occ_coeff*occ_coeff',
%   and the ECP-consistent nuclear charges (see READ_NUCLEAR_CHARGES).
    required = {'filename', 'mol', 'Nbasis', 'Nbasis_indep', 'alpha_MO_coeff', 'xyz_bohr', 'Nalpha', 'Nbeta', 'Nelec'};
    for k = 1:numel(required)
        if ~isfield(data, required{k})
            error('G_draw_esp_surface: %s must be the struct returned by G09_fchk_read/G16_fchk_read (missing "%s" field).', argname, required{k});
        end
    end
    if ~isfile(data.filename)
        error('G_draw_esp_surface: %s.filename (%s) not found -- this function re-reads raw sections directly from the .fchk file.', argname, data.filename);
    end
    if data.Nalpha ~= data.Nbeta
        error('G_draw_esp_surface: %s is an open-shell calculation (Nalpha=%d, Nbeta=%d) -- only closed-shell (Nalpha==Nbeta) is currently supported.', argname, data.Nalpha, data.Nbeta);
    end

    aobasis = g_read_aobasis_from_fchk(data.filename);
    if isempty(aobasis.shell_types)
        error('G_draw_esp_surface: no basis-set sections found in %s (%s) -- is this a valid .fchk file?', data.filename, argname);
    end
    g_check_supported_shells(aobasis.shell_types);

    alpha_MO  = reshape(data.alpha_MO_coeff, data.Nbasis, data.Nbasis_indep);
    occ_coeff = alpha_MO(:, 1:data.Nalpha);
    P_density = 2 * (occ_coeff * occ_coeff');

    nuclear_charges = read_nuclear_charges(data.filename, data.mol.Natoms);
end

function V_esp = compute_esp_at_points(aobasis, occ_coeff, P_density, nuclear_charges, ...
    xyz_bohr, Natoms, Nelec, verts_bohr, esp_method, esp_spacing, esp_pad, dist_floor, verbose, tag)
%COMPUTE_ESP_AT_POINTS  V_nuc - V_elec at VERTS_BOHR for one calculation
%   (either DATA or a 'CompareTo' struct, identified by TAG for verbose
%   messages/error rethrows). Nuclear term is exact/analytic; electronic
%   term is G_EVAL_ESP_ANALYTIC by default, falling back to the numerical
%   Coulomb-grid-sum method -- see the "Method"/"Accuracy" notes in this
%   file's help text, which this function implements unchanged from the
%   single-calculation case, just parameterized so 'CompareTo' mode can
%   call it twice (once per calculation, at the SAME vertices).
    Nv = size(verts_bohr, 1);

    V_nuc = zeros(Nv, 1);
    for A = 1:Natoms
        dd = sqrt(sum((verts_bohr - xyz_bohr(A,:)).^2, 2));
        V_nuc = V_nuc + nuclear_charges(A) ./ max(dd, 1e-6);
    end

    used_analytic = false;
    if ~strcmp(esp_method, 'numeric')
        try
            if verbose
                fprintf('G_draw_esp_surface (%s): evaluating electronic ESP term analytically (McMurchie-Davidson)...\n', tag);
            end
            tic;
            V_elec = g_eval_esp_analytic(aobasis, P_density, verts_bohr, verbose);
            if verbose
                fprintf('  Analytic electronic-term evaluation: %.1f s\n', toc);
            end
            used_analytic = true;
        catch ME
            if strcmp(esp_method, 'analytic') || ~strcmp(ME.identifier, 'g_eval_esp_analytic:unsupportedShell')
                rethrow(ME);
            end
            if verbose
                fprintf('G_draw_esp_surface (%s): basis has shells unsupported by the analytic ESP method (%s) -- falling back to the numerical Coulomb-grid-sum method.\n', tag, ME.message);
            end
        end
    end

    if ~used_analytic
        lo_e = min(xyz_bohr, [], 1) - esp_pad;
        hi_e = max(xyz_bohr, [], 1) + esp_pad;
        gxe = lo_e(1):esp_spacing:hi_e(1);
        gye = lo_e(2):esp_spacing:hi_e(2);
        gze = lo_e(3):esp_spacing:hi_e(3);
        [Xe, Ye, Ze] = meshgrid(gxe, gye, gze);
        esp_grid_pts = [Xe(:), Ye(:), Ze(:)];
        dVe = esp_spacing^3;

        if verbose
            fprintf('G_draw_esp_surface (%s): Coulomb-sum grid %d x %d x %d = %d points (spacing=%.2f Bohr, padding=%.1f Bohr)\n', ...
                tag, numel(gxe), numel(gye), numel(gze), size(esp_grid_pts,1), esp_spacing, esp_pad);
        end

        tic;
        rho_e = g_eval_density_on_grid(aobasis, occ_coeff, esp_grid_pts);
        if verbose
            fprintf('  Coulomb-grid density evaluation: %.1f s\n', toc);
        end

        % See the "Accuracy" note in this file's help text: a uniform grid
        % systematically mis-estimates the raw density at each nuclear
        % cusp, so the grid is rescaled to the exactly-known electron
        % count before the Coulomb sum.
        captured_electrons = sum(rho_e) * dVe;
        scale_factor = Nelec / captured_electrons;
        if verbose
            fprintf('  Electrons captured on Coulomb grid: %.2f / %d -- rescaling density by %.4f\n', ...
                captured_electrons, Nelec, scale_factor);
        end
        rho_e = rho_e * scale_factor;

        V_elec = zeros(Nv, 1);
        batch_size = 100;
        n_batches = ceil(Nv / batch_size);
        tic;
        for b = 1:n_batches
            idx0 = (b-1)*batch_size + 1;
            idx1 = min(b*batch_size, Nv);
            vb = verts_bohr(idx0:idx1, :);   % [Nb x 3]

            dx = esp_grid_pts(:,1) - vb(:,1)';   % [Ngrid x Nb]
            dy = esp_grid_pts(:,2) - vb(:,2)';
            dz = esp_grid_pts(:,3) - vb(:,3)';
            dist = sqrt(dx.^2 + dy.^2 + dz.^2);
            dist = max(dist, dist_floor);

            V_elec(idx0:idx1) = sum((rho_e * dVe) ./ dist, 1)';

            if verbose && (mod(b, max(1,round(n_batches/5))) == 0 || b == n_batches)
                fprintf('  Coulomb sum: batch %d/%d (%.1f s elapsed)\n', b, n_batches, toc);
            end
        end
    end % ~used_analytic

    V_esp = V_nuc - V_elec;
end
