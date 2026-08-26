function rho = g_eval_density_on_grid(aobasis, occ_coeff, grid_pts_bohr, occ_factor)
%G_EVAL_DENSITY_ON_GRID  rho(r) = OCC_FACTOR * sum_i psi_i(r)^2 over the
%   occupied orbitals in OCC_COEFF, using G_EVAL_MO_ON_GRID to evaluate
%   every occupied orbital in a single pass through the basis shells.
%   Shared by G_draw_density_surface and G_draw_esp_surface.
%
%   rho = g_eval_density_on_grid(aobasis, occ_coeff, grid_pts_bohr)
%   rho = g_eval_density_on_grid(aobasis, occ_coeff, grid_pts_bohr, occ_factor)
%
%   occ_coeff: [Nbasis x Nocc] occupied MO coefficient columns (e.g.
%   reshape(data.alpha_MO_coeff, data.Nbasis, data.Nbasis_indep)(:,1:data.Nalpha)).
%   occ_factor: occupation number per orbital column (default: 2, closed-
%   shell double occupancy); pass 1 for a single alpha- or beta-spin
%   channel (e.g. spin density = rho_alpha(occ_factor=1) -
%   rho_beta(occ_factor=1), each spin-orbital holding exactly one
%   electron). Returns rho: [Ngrid x 1], in electrons/Bohr^3.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

if nargin < 4
    occ_factor = 2;
end

psi_occ = g_eval_mo_on_grid(aobasis, occ_coeff, grid_pts_bohr);   % [Ngrid x Nocc]
rho = occ_factor * sum(psi_occ.^2, 2);

end % g_eval_density_on_grid
