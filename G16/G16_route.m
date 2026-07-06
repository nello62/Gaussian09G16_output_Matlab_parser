function route = G16_route(filename)
% G16_ROUTE  Extracts the route section from a Gaussian 16 .out/.log file.
%
%   route = G16_ROUTE(filename)
%
%   OUTPUT  char — full route section string (on a single line)
%
%   Example:
%       r = G16_route('V_E00t.out')
%       % r = '# opt=calcall freq=raman CPHF=Rdfreq b3lyp/6-311g(d,p) nosymm cphf=grid=fine int=grid=ultrafine'
%
%   Note: collects the lines between the two '------' separators that follow
%         the first '#' line (i.e. the standard Gaussian route block).

if ~isfile(filename)
    error('G16_route: file not found: %s', filename);
end

fid  = fopen(filename, 'r');
raw  = fread(fid, '*char')';
fclose(fid);
lines = strsplit(raw, newline);
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
            route_lines{end+1} = ln; %#ok<AGROW>
        else
            % Not a route line (e.g. title before the real separator): reset
            found_first_sep = false;
        end
        continue
    end

    if in_route
        route_lines{end+1} = ln; %#ok<AGROW>
    end
end

if isempty(route_lines)
    error('G16_route: route section not found in %s', filename);
end

% Join all continuation lines into one string
route = strtrim(strjoin(route_lines, ' '));

fprintf('Route: %s\n', route);

end
