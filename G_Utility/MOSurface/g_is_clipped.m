function tf = g_is_clipped(fv, Xa, Ya, Za, tol)
%G_IS_CLIPPED  True if any isosurface vertex lies within TOL of the
%   evaluated grid's outer boundary -- a sign that the true isosurface
%   extends beyond the grid and has been artificially cut off there.
%   Shared by G_draw_mo_surface and G_draw_density_surface.
%
%   tf = g_is_clipped(fv, Xa, Ya, Za, tol)
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

if isempty(fv.vertices)
    tf = false;
    return
end
lo = [min(Xa(:)) min(Ya(:)) min(Za(:))];
hi = [max(Xa(:)) max(Ya(:)) max(Za(:))];
v  = fv.vertices;
tf = any(any(v <= lo + tol | v >= hi - tol));

end % g_is_clipped
