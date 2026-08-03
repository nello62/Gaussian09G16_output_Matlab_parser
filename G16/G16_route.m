function route = G16_route(filename, varargin)
% G16_ROUTE  Extracts the route section from a Gaussian 16 .out/.log file.
%
%   route = G16_ROUTE(filename)
%
%   OUTPUT  char — full route section string (on a single line)
%
%   Optional parameters:
%       'Lines'  - pre-read cell array of file lines, to skip re-reading
%                  the file when it has already been read elsewhere (e.g.
%                  G16_READ_ALL). Default {} (read the file normally).
%
%   Example:
%       r = G16_route('V_E00t.out')
%       % r = '# opt=calcall freq=raman CPHF=Rdfreq b3lyp/6-311g(d,p) nosymm cphf=grid=fine int=grid=ultrafine'
%
%   Note: collects the lines between the two '------' separators that follow
%         the first '#' line (i.e. the standard Gaussian route block).
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

p = inputParser;
addRequired(p,  'filename', @ischar);
addParameter(p, 'Lines',    {}, @iscell);
parse(p, filename, varargin{:});

lines = p.Results.Lines;
if isempty(lines)
    if ~isfile(filename)
        error('G16_route: file not found: %s', filename);
    end
    fid  = fopen(filename, 'r');
    raw  = fread(fid, '*char')';
    fclose(fid);
    lines = strsplit(raw, newline);
    G16_check_gaussian_match(lines, filename);
end
N = numel(lines);

% The route section in G16 is delimited by two '----' lines
% The first '#' line is found between these two separators.
% Schema:
%   --------------...
%   # opt=calcall freq=raman ...
%   continued route...
%   --------------...

route_lines = {};
in_route    = false;
found_first_sep = false;
sep_count   = 0;

for k = 1:N
    ln = strtrim(lines{k});

    % Separator: a line consisting entirely of dashes (>= 20 chars)
    is_sep = ~isempty(ln) && all(ln == '-') && numel(ln) >= 20;

    if is_sep
        sep_count = sep_count + 1;
        if in_route
            % Second separator reached: end of route section
            break
        end
        % First separator found: the following lines should be the route
        found_first_sep = true;
        continue
    end

    if found_first_sep && ~in_route
        % First non-separator line after the first sep: must start with #
        if ~isempty(regexp(ln, '^#', 'once'))
            in_route = true;
            % Strip only the leading indent space; keep any trailing
            % whitespace exactly as printed (see join comment below).
            route_lines{end+1} = regexprep(lines{k}, '^\s+', ''); %#ok<AGROW>
        else
            % Not a route line (e.g. title before the real separator): reset
            found_first_sep = false;
        end
        continue
    end

    if in_route
        route_lines{end+1} = regexprep(lines{k}, '^\s+', ''); %#ok<AGROW>
    end
end

if isempty(route_lines)
    error('G16_route: route section not found in %s', filename);
end

% Gaussian wraps the route echo at a fixed column with no regard for word
% boundaries, so a keyword can be split mid-word across two lines (e.g.
% "nosym" / "m cphf=..." for "nosymm cphf=..."). Joining with an inserted
% space (the previous behaviour) corrupted every such keyword. Since
% Gaussian never leaves a genuine trailing space before the forced wrap,
% concatenating the lines directly (no separator) reconstructs the
% original text correctly in both cases: a mid-word split rejoins
% cleanly, and a wrap that happens to fall between two words still has
% its separating space, because that space is part of the line content
% itself, not the join. A final whitespace-collapse tidies up the rare
% double space this can otherwise leave (e.g. if a wrap coincides with a
% word boundary and both leading/trailing spaces would survive).
route = regexprep(strtrim(strjoin(route_lines, '')), '\s+', ' ');

fprintf('Route: %s\n', route);

end
