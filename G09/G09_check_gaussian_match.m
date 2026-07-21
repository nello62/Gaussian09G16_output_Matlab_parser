function ok = G09_check_gaussian_match(lines, filename)
% G09_CHECK_GAUSSIAN_MATCH  Warns if FILENAME looks like a different
%                           Gaussian major version than this toolbox expects.
%
%   ok = G09_CHECK_GAUSSIAN_MATCH(lines, filename)
%
%   The G09_*.m toolbox functions expect Gaussian 09 output; Gaussian 09
%   and Gaussian 16 differ in output formatting, so using the wrong
%   toolbox on a file can silently misparse it. This scans the
%   already-read LINES (no extra disk I/O) for the "Gaussian NN,
%   Revision X.YY," citation line printed near the top of every output
%   file; if NN is found and is not 9, a non-blocking warning is issued
%   suggesting the G16_*.m toolbox instead. Silently does nothing
%   (returns true) if the citation line is not found near the top of the
%   file.
%
%   OUTPUT  OK - true if no mismatch was detected (or the version is
%                unknown), false if a mismatch warning was issued.
%
%   See also G09_GAUSSIAN_VERSION, G16_CHECK_GAUSSIAN_MATCH.

ok = true;
nscan = min(numel(lines), 60);
for k = 1:nscan
    tok = regexp(lines{k}, 'Gaussian\s+(\d+),\s*Revision\s+([\w.]+)', 'tokens', 'once');
    if ~isempty(tok)
        major = str2double(tok{1});
        if major ~= 9
            warning('G09_check_gaussian_match:versionMismatch', ...
                ['%s looks like a Gaussian %s output file, not Gaussian 09.\n' ...
                 'Consider using the G16_*.m toolbox functions instead ' ...
                 '(G09 and G16 output formats differ).'], ...
                filename, tok{1});
            ok = false;
        end
        return
    end
end

end % G09_check_gaussian_match
