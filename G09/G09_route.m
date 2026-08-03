function route = G09_route(filename)
% G09_ROUTE  Extracts the route section from a Gaussian 09 output file.
%
%   route = G09_ROUTE(filename)
%
%   Returns the complete route section as a single string.
%   The format is identical to G16 (--- separator pair around '#' lines).
%
%   Example:
%       r = G09_route('indaco.log')
%       % '# opt freq=raman b3lyp/6-311++g(d,p) nosymm geom=connectivity field=x+50'
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

if ~isfile(filename)
    error('G09_route: file not found: %s', filename);
end

lines = G09_read_lines(filename);
N     = numel(lines);

% Find separator that is followed within 3 lines by a '#' line
sep_idx = find(~cellfun(@isempty, regexp(lines, '^\s*-{10,}\s*$')));

route_sep_start = [];
for si = 1:numel(sep_idx)
    k = sep_idx(si);
    for j = k+1 : min(k+3, N)
        if ~isempty(regexp(strtrim(lines{j}), '^#', 'once'))
            route_sep_start = k;
            break
        end
    end
    if ~isempty(route_sep_start), break; end
end

if isempty(route_sep_start)
    error('G09_route: route section not found in %s', filename);
end

route_sep_end = [];
for si = 1:numel(sep_idx)
    if sep_idx(si) > route_sep_start
        route_sep_end = sep_idx(si);
        break
    end
end

if isempty(route_sep_end)
    error('G09_route: closing separator for route not found in %s', filename);
end

route_lines = {};
for k = route_sep_start+1 : route_sep_end-1
    ln = strtrim(lines{k});
    if ~isempty(ln)
        % Strip only the leading indent space; keep any trailing
        % whitespace exactly as printed (see join comment below).
        route_lines{end+1} = regexprep(lines{k}, '^\s+', ''); %#ok<AGROW>
    end
end

% Gaussian wraps the route echo at a fixed column with no regard for word
% boundaries, so a keyword can be split mid-word across two lines (e.g.
% "nosym" / "m geom=..." for "nosymm geom=..."). Joining with an inserted
% space (the previous behaviour) corrupted every such keyword. Since
% Gaussian never leaves a genuine trailing space before the forced wrap,
% concatenating the lines directly (no separator) reconstructs the
% original text correctly in both cases: a mid-word split rejoins
% cleanly, and a wrap that happens to fall between two words still has
% its separating space, because that space is part of the line content
% itself, not the join. A final whitespace-collapse tidies up the rare
% double space this can otherwise leave.
route = regexprep(strtrim(strjoin(route_lines, '')), '\s+', ' ');
fprintf('Route: %s\n', route);

end  % G09_route
