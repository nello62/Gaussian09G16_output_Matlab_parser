function ginp = G_read_input(filename)
% G_READ_INPUT  Reads a Gaussian input file (.gjf/.com/.in) and extracts
%               link0 commands, route section, title, charge/multiplicity,
%               and the starting Cartesian geometry.
%
%   ginp = G_READ_INPUT(filename)
%
%   The Gaussian input file format does not differ between Gaussian 09 and
%   Gaussian 16, so this is a single shared function rather than separate
%   G09_/G16_ versions -- an identical copy is kept in both the G09/ and
%   G16/ folders, so it is available regardless of which toolbox you have
%   on the MATLAB path.
%
%   OUTPUT  struct ginp with fields:
%       .symbols    {Natoms x 1 cell}   atomic symbols
%       .xyz        [Natoms x 3 double] starting Cartesian coordinates (Angstrom)
%       .Natoms     int
%       .charge     int      total molecular charge
%       .mult       int      spin multiplicity
%       .title      char     title/comment line(s), joined with a space
%       .route      char     full route section, single line
%       .method     char     method parsed from the route (e.g. 'B3LYP'), '' if not found
%       .basis      char     basis set parsed from the route (e.g. '6-311+G(d,p)'), '' if not found
%       .chk        char     %chk link0 value, '' if absent
%       .mem        char     %mem link0 value, '' if absent
%       .nproc      char     %nprocshared/%nproc link0 value, '' if absent
%       .filename   char     source file path
%
%   ginp has the same .symbols/.xyz/.Natoms/.filename fields as the struct
%   returned by G09_STRUCTURE/G16_STRUCTURE, so it can be passed directly
%   to G09_DRAW_MOLECULE/G16_DRAW_MOLECULE or G09_GET_BOND_LENGTH/
%   G16_GET_BOND_LENGTH without any conversion.
%
%   Example:
%       ginp = G_read_input('molecule.gjf');
%       G16_draw_molecule(ginp, 'ShowAxes', true);
%       fprintf('%s/%s, charge %d, mult %d\n', ginp.method, ginp.basis, ...
%           ginp.charge, ginp.mult);
%
%   Limitation: only Cartesian-coordinate geometry blocks are supported
%   (not Z-matrix input).
%
%   See also G09_STRUCTURE, G16_STRUCTURE, G09_DRAW_MOLECULE,
%            G16_DRAW_MOLECULE, G09_RESTART.

if ~isfile(filename)
    error('G_read_input: file not found: %s', filename);
end

fid = fopen(filename, 'r');
raw = fread(fid, '*char')';
fclose(fid);
raw = strrep(raw, sprintf('\r\n'), newline);
raw = strrep(raw, sprintf('\r'), newline);
lines = strsplit(raw, newline, 'CollapseDelimiters', false);

% Drop pure comment lines ('!' as the first non-blank character), which
% Gaussian allows anywhere in an input file and which carry no structural
% meaning for this parser.
is_comment = cellfun(@(l) ~isempty(regexp(strtrim(l), '^!', 'once')), lines);
lines = lines(~is_comment);
N = numel(lines);

is_separator = @(s) ~isempty(regexp(s, '^-{4,}$', 'once'));

% -------------------------------------------------------------------------
% Link0 (%) commands
% -------------------------------------------------------------------------
chk = ''; mem = ''; nproc = '';
i = 1;
while i <= N
    ln = strtrim(lines{i});
    if isempty(ln)
        i = i + 1;
        continue
    end
    if ln(1) ~= '%'
        break
    end
    tok = regexp(ln, '^%(\w+)\s*=\s*(.+)$', 'tokens', 'once');
    if ~isempty(tok)
        key = lower(tok{1});
        val = strtrim(tok{2});
        switch key
            case 'chk',                      chk   = val;
            case 'mem',                       mem   = val;
            case {'nprocshared', 'nproc'},    nproc = val;
        end
    end
    i = i + 1;
end

% Skip a decorative separator line before the route, if present
while i <= N && is_separator(strtrim(lines{i}))
    i = i + 1;
end
while i <= N && isempty(strtrim(lines{i}))
    i = i + 1;
end

% -------------------------------------------------------------------------
% Route section: from the first '#' line to the next blank line
% -------------------------------------------------------------------------
route_lines = {};
first_trimmed = '';
if i <= N
    first_trimmed = strtrim(lines{i});
end
if ~isempty(first_trimmed) && first_trimmed(1) == '#'
    while i <= N
        tln = strtrim(lines{i});
        if isempty(tln) || is_separator(tln)
            % A blank line or a closing separator both end the route:
            % some generated inputs (e.g. G09_restart) bracket the route
            % between two separator lines with no blank line before the
            % title that follows, so a separator must terminate the
            % route immediately rather than being skipped over.
            i = i + 1;
            break
        end
        route_lines{end+1} = tln; %#ok<AGROW>
        i = i + 1;
    end
end
route = strtrim(strjoin(route_lines, ' '));

method = ''; basis = '';
rtoks = strsplit(route);
mb_idx = find(contains(rtoks, '/'), 1);
if ~isempty(mb_idx)
    parts = strsplit(rtoks{mb_idx}, '/');
    if numel(parts) >= 2
        method = parts{1};
        basis  = parts{2};
    end
end

% Skip a decorative separator line after the route, and blank lines
while i <= N && (isempty(strtrim(lines{i})) || is_separator(strtrim(lines{i})))
    i = i + 1;
end

% -------------------------------------------------------------------------
% Title (one or more non-blank lines, until the next blank line)
% -------------------------------------------------------------------------
title_lines = {};
while i <= N && ~isempty(strtrim(lines{i}))
    title_lines{end+1} = strtrim(lines{i}); %#ok<AGROW>
    i = i + 1;
end
titlestr = strjoin(title_lines, ' ');
i = i + 1;   % skip the blank line after the title

% -------------------------------------------------------------------------
% Charge / multiplicity line
% -------------------------------------------------------------------------
charge = NaN; mult = NaN;
while i <= N && isempty(strtrim(lines{i}))
    i = i + 1;
end
if i <= N
    tok = regexp(strtrim(lines{i}), '^(-?\d+)\s+(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        charge = str2double(tok{1});
        mult   = str2double(tok{2});
    end
    i = i + 1;
end

% -------------------------------------------------------------------------
% Geometry: "Symbol  X  Y  Z" Cartesian lines, until a blank line or EOF
% -------------------------------------------------------------------------
symbols = {};
xyz = zeros(0, 3);
while i <= N
    ln = strtrim(lines{i});
    if isempty(ln)
        break
    end
    parts = strsplit(ln);
    if numel(parts) >= 4
        x = str2double(parts{end-2});
        y = str2double(parts{end-1});
        z = str2double(parts{end});
        if ~any(isnan([x y z]))
            symtok = regexp(parts{1}, '^([A-Za-z]{1,2})', 'tokens', 'once');
            if ~isempty(symtok)
                sym = symtok{1};
                if numel(sym) > 1
                    sym = [upper(sym(1)), lower(sym(2))];
                else
                    sym = upper(sym);
                end
                symbols{end+1, 1} = sym;    %#ok<AGROW>
                xyz(end+1, :) = [x, y, z];  %#ok<AGROW>
            end
        end
    end
    i = i + 1;
end

if isempty(symbols)
    error('G_read_input:noGeometry', ...
        'No Cartesian geometry block found in %s (Z-matrix input is not supported).', filename);
end

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
ginp.symbols  = symbols;
ginp.xyz      = xyz;
ginp.Natoms   = numel(symbols);
ginp.charge   = charge;
ginp.mult     = mult;
ginp.title    = titlestr;
ginp.route    = route;
ginp.method   = method;
ginp.basis    = basis;
ginp.chk      = chk;
ginp.mem      = mem;
ginp.nproc    = nproc;
ginp.filename = filename;

fprintf('\n── G_read_input: %s ──\n', filename);
fprintf('  Route  : %s\n', route);
if ~isempty(method)
    fprintf('  Method/Basis : %s / %s\n', method, basis);
end
fprintf('  Charge = %d   Multiplicity = %d\n', charge, mult);
fprintf('  %d atoms\n\n', ginp.Natoms);

end % G_read_input
