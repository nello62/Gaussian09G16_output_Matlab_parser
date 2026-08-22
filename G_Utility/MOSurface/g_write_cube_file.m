function g_write_cube_file(filename, title_line, comment_line, atomic_numbers, xyz_bohr, gx, gy, gz, V, orbital_indices)
%G_WRITE_CUBE_FILE  Writes a scalar volumetric grid to a Gaussian-format
%   .cube file (the same format read/validated earlier against real
%   cubegen output -- see G_EVAL_ESP_ANALYTIC's header and
%   Theory_Vibrations_Polarizability.tex Part III).
%
%   g_write_cube_file(filename, title_line, comment_line, ...
%       atomic_numbers, xyz_bohr, gx, gy, gz, V)
%   g_write_cube_file(..., orbital_indices)   % MO cube: negative Natoms
%                                              % + the extra orbital line
%
%   atomic_numbers : [Natoms x 1], true atomic numbers (NOT the
%                    ECP-reduced "Nuclear charges" -- Gaussian's own
%                    cubegen writes the real atomic number here, only
%                    ESP/density VALUES depend on the ECP-reduced charge)
%   xyz_bohr       : [Natoms x 3]
%   gx, gy, gz     : 1D axis vectors (Bohr) as built by
%                    LO(k):SPACING:HI(k) in the calling G_DRAW_*_SURFACE
%                    function -- assumed uniformly spaced, giving an
%                    axis-aligned (diagonal-step) cube grid
%   V              : scalar values in MESHGRID convention, i.e. the SAME
%                    array produced by
%                        [X,Y,Z] = meshgrid(gx,gy,gz);
%                        V = reshape(f(X(:),Y(:),Z(:)), size(X));
%                    as already used for ISOSURFACE/ISONORMALS by every
%                    G_DRAW_*_SURFACE caller -- size [numel(gy),
%                    numel(gx), numel(gz)], i.e. V(iy,ix,iz). This
%                    function handles the meshgrid-to-cube-file-order
%                    permutation internally (verified against an
%                    independent nested-loop reference; see the "Value
%                    ordering" note below), so callers never need to
%                    permute V themselves.
%   orbital_indices- (optional) vector of 1-based MO indices. If given
%                    (and nonempty), Natoms is written as NEGATIVE and an
%                    extra line "Norbitals idx1 idx2 ..." follows the
%                    atom list, per the Gaussian MO-cube convention
%                    (confirmed against a real cubegen-generated MO cube
%                    earlier this session).
%
%   Value ordering: a Gaussian cube file's flat data section is X
%   slowest-varying, Z fastest-varying. Given V in meshgrid convention
%   (V(iy,ix,iz)), the correct flat order is
%   reshape(permute(V,[3 1 2]), [], 1) -- verified by direct construction
%   against an explicit triple nested-loop (ix outer, iy, iz inner)
%   reference on a small test grid, exact match.
%
%   Line wrapping: values are written 6 per line, but -- confirmed
%   against real cubegen output, NOT just assumed -- the line is also
%   force-broken at the end of EVERY z-run (each fixed (ix,iy) pair's
%   Nz values), even if that leaves a short last line, with the next
%   z-run always starting its own fresh line (i.e. z-runs never share a
%   line). Wrapping only every 6 across the whole flat stream, ignoring
%   this z-run boundary, is silently accepted by this toolbox's own
%   token-based cube reader (fscanf ignores line breaks entirely) but
%   was found -- the hard way, from a user report -- to badly confuse
%   GaussView's own cube import: it would render the surface once and
%   then crash, rather than failing to parse outright.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

if nargin < 10
    orbital_indices = [];
end

Natoms = numel(atomic_numbers);
Nx = numel(gx); Ny = numel(gy); Nz = numel(gz);
if ~isequal(size(V), [Ny, Nx, Nz])
    error('g_write_cube_file:sizeMismatch', ...
        'V has size %s, expected [%d %d %d] = [numel(gy) numel(gx) numel(gz)] (meshgrid convention).', ...
        mat2str(size(V)), Ny, Nx, Nz);
end

dx = local_step(gx); dy = local_step(gy); dz = local_step(gz);
origin_bohr = [gx(1), gy(1), gz(1)];

fid = fopen(filename, 'w');
if fid < 0
    error('g_write_cube_file:cannotOpen', 'Could not open %s for writing.', filename);
end
c = onCleanup(@() fclose(fid));

fprintf(fid, '%s\n', title_line);
fprintf(fid, '%s\n', comment_line);

is_mo_cube = ~isempty(orbital_indices);
natoms_field = Natoms;
if is_mo_cube
    natoms_field = -Natoms;
end
fprintf(fid, '%5d%12.6f%12.6f%12.6f\n', natoms_field, origin_bohr(1), origin_bohr(2), origin_bohr(3));
fprintf(fid, '%5d%12.6f%12.6f%12.6f\n', Nx, dx, 0, 0);
fprintf(fid, '%5d%12.6f%12.6f%12.6f\n', Ny, 0, dy, 0);
fprintf(fid, '%5d%12.6f%12.6f%12.6f\n', Nz, 0, 0, dz);

for A = 1:Natoms
    fprintf(fid, '%5d%12.6f%12.6f%12.6f%12.6f\n', atomic_numbers(A), atomic_numbers(A), ...
        xyz_bohr(A,1), xyz_bohr(A,2), xyz_bohr(A,3));
end

if is_mo_cube
    fprintf(fid, '%5d', numel(orbital_indices));
    fprintf(fid, '%5d', orbital_indices(:)');
    fprintf(fid, '\n');
end

% Flatten to Gaussian cube order (X slowest, Z fastest) -- see the
% "Value ordering" note above. permute(V,[3 1 2]) has size [Nz,Ny,Nx]
% (dim1=Z fastest); reshaping to [Nz x (Ny*Nx)] gives one full z-run per
% COLUMN, in the correct overall (X slowest, Y, Z fastest) column order.
zruns = reshape(permute(V, [3 1 2]), Nz, []);

% A real cube file force-breaks the line at the END OF EVERY Z-RUN, not
% just every 6 values continuously across the whole flat stream --
% confirmed against real cubegen output (e.g. a Z-dimension of 104 wraps
% as seventeen 6-value lines then one 2-value line, PER z-run, with the
% next run always starting its own fresh line). Wrapping only every 6
% across the whole file (ignoring the z-run boundary) is silently
% accepted by this toolbox's own token-based reader (fscanf ignores line
% breaks) but was found to confuse GaussView's cube import badly enough
% to crash it after an initially-successful-looking surface render.
for r = 1:size(zruns, 2)
    col = zruns(:, r);
    for k0 = 1:6:Nz
        k1 = min(k0+5, Nz);
        fprintf(fid, '%13.5E', col(k0:k1));
        fprintf(fid, '\n');
    end
end

end % g_write_cube_file

function d = local_step(g)
%LOCAL_STEP  Axis spacing from a 1D grid vector; arbitrary (unused by any
%   real geometry) for a degenerate single-point axis.
    if numel(g) > 1
        d = g(2) - g(1);
    else
        d = 1;
    end
end
