function g_check_supported_shells(shell_types)
%G_CHECK_SUPPORTED_SHELLS  Raises a clear error (no partial rendering) if
%   the basis contains any shell type beyond S/P/SP/D/F/G (pure or
%   Cartesian). Shared by G_draw_mo_surface and G_draw_density_surface.
%
%   g_check_supported_shells(shell_types)
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

supported = [0, 1, -1, -2, 2, -3, 3, -4, 4];
bad = find(~ismember(shell_types, supported));
if isempty(bad)
    return
end
bad_types = unique(shell_types(bad));
names = cell(1, numel(bad_types));
for i = 1:numel(bad_types)
    names{i} = sprintf('%d', bad_types(i));
end
error('g_check_supported_shells:unsupportedShell', ...
    ['This basis set contains shell type(s) [%s] not supported by this function ' ...
     '(only S, P, SP, D, F, and G -- pure or Cartesian -- shells are implemented). ' ...
     'Rendering would be silently wrong for any basis function belonging to these ' ...
     'shells, so no surface is produced. A common cause is an h-polarized or larger ' ...
     'basis set (rare) or a non-standard shell-type code.'], ...
    strjoin(names, ', '));

end % g_check_supported_shells
