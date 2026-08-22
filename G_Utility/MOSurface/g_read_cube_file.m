function [atomic_numbers, xyz_bohr, gx, gy, gz, V, orbital_indices, title_line, comment_line] = g_read_cube_file(filename)
%G_READ_CUBE_FILE  Reads a Gaussian-format .cube file (the inverse of
%   G_WRITE_CUBE_FILE, and compatible with real cubegen/GaussView output,
%   not just files this toolbox wrote itself).
%
%   [atomic_numbers, xyz_bohr, gx, gy, gz, V, orbital_indices, ...
%       title_line, comment_line] = g_read_cube_file(filename)
%
%   atomic_numbers : [Natoms x 1]
%   xyz_bohr       : [Natoms x 3]
%   gx, gy, gz     : 1D axis vectors (Bohr), built from the cube's origin
%                    and (diagonal) step vectors -- a non-axis-aligned
%                    (non-diagonal step matrix) cube raises a clear
%                    error, since this toolbox's grid-based tools
%                    (isosurface/meshgrid) all assume an axis-aligned
%                    grid
%   V              : scalar values in MESHGRID convention, i.e. the SAME
%                    convention G_WRITE_CUBE_FILE expects and every
%                    G_DRAW_*_SURFACE caller already uses internally --
%                    size [numel(gy), numel(gx), numel(gz)], V(iy,ix,iz).
%                    Reconstructed from the file's flat (X slowest, Z
%                    fastest) order via the same reshape/permute recipe
%                    validated earlier against real cubegen output
%                    (relative difference ~1e-7, text-precision-limited)
%   orbital_indices: [] for a scalar-field cube; the listed MO index/
%                    indices for an MO cube (negative Natoms header field)
%   title_line, comment_line: the file's first two free-text lines, used
%                    by G_DRAW_CUBE_SURFACE to guess the field type
%                    (density/MO/ESP) for a sensible default 'IsoValue'
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

fid = fopen(filename, 'r');
if fid < 0
    error('g_read_cube_file:cannotOpen', 'Could not open %s for reading.', filename);
end
c = onCleanup(@() fclose(fid));

title_line   = fgetl(fid);
comment_line = fgetl(fid);

l3 = sscanf(fgetl(fid), '%f')';
Natoms_field = round(l3(1));
origin_bohr  = l3(2:4);
is_mo_cube   = Natoms_field < 0;
Natoms       = abs(Natoms_field);

N = zeros(1,3); step = zeros(3,3);
for k = 1:3
    lk = sscanf(fgetl(fid), '%f')';
    N(k)      = round(lk(1));
    step(k,:) = lk(2:4);
end
if any(abs(step - diag(diag(step))) > 1e-9, 'all')
    error('g_read_cube_file:nonAxisAligned', ...
        '%s has a non-axis-aligned grid (off-diagonal step vectors) -- not supported by this toolbox''s grid-based tools.', filename);
end

atomic_numbers = zeros(Natoms,1);
xyz_bohr = zeros(Natoms,3);
for a = 1:Natoms
    la = sscanf(fgetl(fid), '%f')';
    atomic_numbers(a) = round(la(1));
    xyz_bohr(a,:) = la(3:5);
end

orbital_indices = [];
if is_mo_cube
    ln = sscanf(fgetl(fid), '%f')';
    n_orb = round(ln(1));
    orbital_indices = round(ln(2:1+n_orb));
end

vals = fscanf(fid, '%f');
Nx = N(1); Ny = N(2); Nz = N(3);
if numel(vals) ~= Nx*Ny*Nz
    error('g_read_cube_file:sizeMismatch', ...
        '%s: read %d data values, expected %d (%d x %d x %d).', ...
        filename, numel(vals), Nx*Ny*Nz, Nx, Ny, Nz);
end

% Flat order is X slowest, Z fastest -- reshape as [Nz,Ny,Nx] (column-
% major, Z fastest) then permute to the natural (ix,iy,iz) order, exactly
% the recipe validated against real cubegen output earlier this session.
Vxyz = permute(reshape(vals, Nz, Ny, Nx), [3 2 1]);   % Vxyz(ix,iy,iz)
V = permute(Vxyz, [2 1 3]);                            % meshgrid convention V(iy,ix,iz)

gx = origin_bohr(1) + (0:Nx-1) * step(1,1);
gy = origin_bohr(2) + (0:Ny-1) * step(2,2);
gz = origin_bohr(3) + (0:Nz-1) * step(3,3);

end % g_read_cube_file
