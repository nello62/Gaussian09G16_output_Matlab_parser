function params = g_cube_grid_params(cubefile)
%G_CUBE_GRID_PARAMS  Recovers the grid spacing and padding used to build
%   a .cube file, so a companion cube (e.g. for G_DRAW_CUBE_SURFACE's
%   'ColorBy', which requires two cubes to share an identical grid) can
%   be regenerated with a matching 'GridSpacing'/'Padding' (or
%   'CubeSpacing'/'CubePadding' for G_DRAW_ESP_SURFACE's own SaveCube
%   grid) without having to note them down manually.
%
%   params = g_cube_grid_params(cubefile)
%
%   Works on any axis-aligned, uniformly-spaced .cube file containing
%   atoms -- one written by this toolbox's G_WRITE_CUBE_FILE (via any
%   G_DRAW_*_SURFACE 'SaveCube' option), or a real cubegen/GaussView
%   cube.
%
%   Output: struct with fields
%       GridSpacing     - scalar, Bohr (the cube's step size; errors if
%                         dx/dy/dz differ by more than 1e-6 Bohr, i.e.
%                         the cube is not isotropic -- this toolbox
%                         always writes isotropic cubes, so a mismatch
%                         means a non-toolbox cube)
%       Padding         - scalar, Bohr: the mean of the six per-face
%                         paddings below -- the value to pass back as
%                         'Padding'/'CubePadding' for a matching cube
%       PaddingPerFace  - [2x3], Bohr: [low;high] x [x,y,z], the raw
%                         per-face (grid boundary minus atom bounding
%                         box) distances Padding was averaged from. A
%                         warning is issued if these differ by more than
%                         one GridSpacing, since this toolbox's own
%                         grids always use the SAME padding on every
%                         side -- a wide spread here means the cube was
%                         not built that way, and 'Padding' alone will
%                         not exactly reproduce this grid.
%       Npoints         - [Nx Ny Nz]
%       Origin_bohr     - [1x3]
%
%   CAUTION -- direction matters: ESP is far more expensive to evaluate
%   per grid point than density or an MO (McMurchie-Davidson integrals,
%   or a Coulomb-grid-sum fallback for basis sets with shells the
%   analytic method does not support), and G_DRAW_ESP_SURFACE's own
%   'SaveCube' path evaluates ESP at every point of the saved grid, not
%   just at isosurface vertices, capped at 2e6 points so it errors out
%   instead of silently running for hours. A grid recovered from a
%   DENSITY cube (typically 0.10-0.15 Bohr spacing, fine enough for a
%   nice-looking isosurface) is very often too fine to reuse for an ESP
%   'SaveCube' this way. Prefer the opposite order: pick a shared
%   GridSpacing/Padding the ESP evaluation can afford (e.g. around
%   G_draw_esp_surface's own defaults, CubeSpacing=0.30/CubePadding=4.0)
%   for the ESP cube first, then match the (much cheaper to re-evaluate)
%   density cube to IT:
%
%       G_draw_esp_surface(data, 'SaveCube', 'esp.cube', ...
%           'CubeSpacing', 0.30, 'CubePadding', 4.0);
%       p = g_cube_grid_params('esp.cube');
%       G_draw_density_surface(data, 'SaveCube', 'density.cube', ...
%           'GridSpacing', p.GridSpacing, 'Padding', p.Padding);
%       G_draw_cube_surface('density.cube', 'ColorBy', 'esp.cube');
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

[~, xyz_bohr, gx, gy, gz] = g_read_cube_file(cubefile);

dx = gx(2) - gx(1);
dy = gy(2) - gy(1);
dz = gz(2) - gz(1);
if max(abs([dx dy dz] - dx)) > 1e-6
    error('g_cube_grid_params:anisotropic', ...
        '%s has different spacing along x/y/z (%.6g/%.6g/%.6g Bohr) -- not an isotropic grid as written by this toolbox.', ...
        cubefile, dx, dy, dz);
end
spacing = dx;

if isempty(xyz_bohr)
    error('g_cube_grid_params:noAtoms', ...
        '%s has no atoms -- cannot recover Padding (only GridSpacing).', cubefile);
end

atom_lo = min(xyz_bohr, [], 1);
atom_hi = max(xyz_bohr, [], 1);
grid_lo = [gx(1), gy(1), gz(1)];
grid_hi = [gx(end), gy(end), gz(end)];

pad_lo  = atom_lo - grid_lo;   % [1x3]
pad_hi  = grid_hi - atom_hi;   % [1x3]
pad_all = [pad_lo; pad_hi];    % [2x3]

padding = mean(pad_all(:));
if max(pad_all(:)) - min(pad_all(:)) > spacing
    warning('g_cube_grid_params:nonUniformPadding', ...
        ['%s: the six face paddings vary by more than one grid step ' ...
         '(%.3g to %.3g Bohr) -- this cube may not have been built with ' ...
         'a single uniform ''Padding'' value; the returned Padding ' ...
         '(%.3g Bohr, the mean) may not exactly reproduce this grid.'], ...
        cubefile, min(pad_all(:)), max(pad_all(:)), padding);
end

params.GridSpacing    = spacing;
params.Padding        = padding;
params.PaddingPerFace = pad_all;
params.Npoints        = [numel(gx), numel(gy), numel(gz)];
params.Origin_bohr    = grid_lo;

end
