function ok = G16_check_gaussian_match(lines, filename)
% G16_CHECK_GAUSSIAN_MATCH  Warns if FILENAME looks like a different
%                           Gaussian major version than this toolbox expects.
%
%   ok = G16_CHECK_GAUSSIAN_MATCH(lines, filename)
%
%   The G16_*.m toolbox functions expect Gaussian 16 output; Gaussian 09
%   and Gaussian 16 differ in output formatting, so using the wrong
%   toolbox on a file can silently misparse it. This scans the
%   already-read LINES (no extra disk I/O) for the "Gaussian NN,
%   Revision X.YY," citation line printed near the top of every output
%   file; if NN is found and is not 16, a non-blocking warning is issued
%   suggesting the G09_*.m toolbox instead. Silently does nothing
%   (returns true) if the citation line is not found near the top of the
%   file.
%
%   OUTPUT  OK - true if no mismatch was detected (or the version is
%                unknown), false if a mismatch warning was issued.
%
%   See also G16_GAUSSIAN_VERSION, G09_CHECK_GAUSSIAN_MATCH.

ok = true;
nscan = min(numel(lines), 60);
for k = 1:nscan
    tok = regexp(lines{k}, 'Gaussian\s+(\d+),\s*Revision\s+([\w.]+)', 'tokens', 'once');
    if ~isempty(tok)
        major = str2double(tok{1});
        if major ~= 16
            warning('G16_check_gaussian_match:versionMismatch', ...
                ['%s looks like a Gaussian %s output file, not Gaussian 16.\n' ...
                 'Consider using the G09_*.m toolbox functions instead ' ...
                 '(G09 and G16 output formats differ).'], ...
                filename, tok{1});
            ok = false;
        end
        return
    end
end

end % G16_check_gaussian_match
