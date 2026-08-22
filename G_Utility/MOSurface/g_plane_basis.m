function [center, u, v, n] = g_plane_basis(xyz_bohr, plane_choice)
%G_PLANE_BASIS  Returns the origin CENTER [1x3], an orthonormal in-plane
%   basis U, V [3x1 each], and the plane normal N [3x1] for a
%   'Mode','contour' cutting plane. Shared by G_draw_mo_surface and
%   G_draw_density_surface.
%
%   [center, u, v, n] = g_plane_basis(xyz_bohr, plane_choice)
%
%   plane_choice: 'auto' -- the best-fit plane through all atoms (PCA --
%   the plane minimizing the atoms' summed squared out-of-plane
%   distance), via the SVD of the centered coordinates: the singular
%   vectors for the two largest singular values span the plane, the
%   third (smallest) is its normal N -- used, together with a plane
%   offset, to probe a plane parallel to (but offset from) the fitted
%   molecular plane, which is essential for a pi-symmetry MO: its
%   amplitude is antisymmetric about (i.e. exactly zero on) the molecular
%   plane itself, so an unoffset auto plane shows nothing for it.
%   'xy'/'xz'/'yz': the corresponding coordinate plane, through the
%   molecule's centroid.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

center = mean(xyz_bohr, 1);
switch plane_choice
    case 'auto'
        if size(xyz_bohr, 1) < 3
            u = [1;0;0]; v = [0;1;0]; n = [0;0;1];   % degenerate (<3 atoms): fall back to xy
        else
            [~, ~, Vsvd] = svd(xyz_bohr - center, 0);
            u = Vsvd(:,1);
            v = Vsvd(:,2);
            n = Vsvd(:,3);
        end
    case 'xy'
        u = [1;0;0]; v = [0;1;0]; n = [0;0;1];
    case 'xz'
        u = [1;0;0]; v = [0;0;1]; n = [0;1;0];
    case 'yz'
        u = [0;1;0]; v = [0;0;1]; n = [1;0;0];
end

end % g_plane_basis
