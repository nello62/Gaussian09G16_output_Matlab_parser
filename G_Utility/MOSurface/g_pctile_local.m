function pv = g_pctile_local(v, p)
%G_PCTILE_LOCAL  P-th percentile of the vector V, via a plain sort (no
%   Statistics and Machine Learning Toolbox dependency, consistent with
%   the rest of this toolbox). Shared by G_draw_density_surface and
%   G_draw_esp_surface.
%
%   pv = g_pctile_local(v, p)
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

v = sort(v(:));
n = numel(v);
if n == 0
    pv = 0;
    return
end
idx = max(1, min(n, round(p/100 * n)));
pv = v(idx);

end % g_pctile_local
