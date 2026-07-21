function mol = G09_structure(filename, varargin)
% G09_STRUCTURE  Extracts the molecular geometry from a Gaussian 09 output file.
%
%   mol = G09_STRUCTURE(filename)
%   mol = G09_STRUCTURE(filename, 'step', 'last')
%   mol = G09_STRUCTURE(filename, 'step', 'first')
%   mol = G09_STRUCTURE(filename, 'step', N)
%
%   Gaussian 09 writes only "Input orientation:" blocks (never Standard
%   orientation), so the 'orientation' parameter of G16_structure is not
%   needed here.
%
%   OUTPUT  struct mol with fields:
%       .symbols     {Natoms x 1 cell}    atomic symbols
%       .xyz         [Natoms x 3 double]  coordinates in Angstrom (X Y Z)
%       .Z           [Natoms x 1 int]     atomic numbers
%       .Natoms      int
%       .step        int                  index of extracted step (1-based)
%       .n_steps     int                  total geometry blocks in file
%       .orientation char                 always 'Input orientation'
%       .filename    char
%
%   Optional parameters also include:
%       'Lines'  - pre-read cell array of file lines (from G09_READ_LINES
%                  or this function's own reader), to skip re-reading the
%                  file when it has already been read elsewhere (e.g.
%                  G09_READ_ALL). Default {} (read the file normally).
%
%   Example:
%       mol = G09_structure('indaco.log');
%       mol = G09_structure('indaco.log', 'step', 1);

% -------------------------------------------------------------------------
% Parse input arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename', @ischar);
addParameter(p, 'step',     'last', @(x) ischar(x) || isnumeric(x));
addParameter(p, 'Lines',    {},     @iscell);
parse(p, filename, varargin{:});
step_req = p.Results.step;

% -------------------------------------------------------------------------
% Read file  (G09 uses CRLF + latin-1 encoding)
% -------------------------------------------------------------------------
lines = p.Results.Lines;
if isempty(lines)
    if ~isfile(filename)
        error('G09_structure: file not found: %s', filename);
    end
    fid  = fopen(filename, 'r', 'n', 'ISO-8859-1');
    raw  = fread(fid, '*char')';
    fclose(fid);

    % Normalise line endings: remove \r
    raw   = strrep(raw, sprintf('\r\n'), newline);
    raw   = strrep(raw, sprintf('\r'),   newline);
    lines = strsplit(raw, newline);

    G09_check_gaussian_match(lines, filename);
end
N = numel(lines);

% -------------------------------------------------------------------------
% Z -> symbol table
% -------------------------------------------------------------------------
Z2sym = Z2symbol_table();

% -------------------------------------------------------------------------
% Find all "Input orientation:" blocks
% -------------------------------------------------------------------------
block_starts = find(~cellfun(@isempty, strfind(lines, 'Input orientation:')));

if isempty(block_starts)
    error('G09_structure: no "Input orientation:" block found in %s', filename);
end

n_blocks = numel(block_starts);

% Resolve step
if ischar(step_req)
    switch lower(step_req)
        case 'last',  step_idx = n_blocks;
        case 'first', step_idx = 1;
        otherwise, error('G09_structure: step must be ''first'', ''last'', or an integer.');
    end
else
    step_idx = round(step_req);
    if step_idx < 1 || step_idx > n_blocks
        error('G09_structure: step %d out of range [1, %d].', step_idx, n_blocks);
    end
end

header_line = block_starts(step_idx);
data_start  = header_line + 5;   % skip: sep + col_header1 + col_header2 + sep

% -------------------------------------------------------------------------
% Parse atom data
% -------------------------------------------------------------------------
symbols = {};
XYZ     = [];
Zvec    = [];

for k = data_start : N
    ln = strtrim(lines{k});
    if isempty(ln), continue; end
    if all(ln == '-'), break; end

    m = regexp(ln, '^\s*\d+\s+(\d+)\s+\d+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', ...
               'tokens', 'once');
    if isempty(m), continue; end

    znum = str2double(m{1});
    x    = str2double(m{2});
    y    = str2double(m{3});
    z    = str2double(m{4});

    symbols{end+1, 1} = Z2sym(int32(znum)); %#ok<AGROW>
    XYZ(end+1, :)     = [x, y, z];          %#ok<AGROW>
    Zvec(end+1, 1)    = znum;               %#ok<AGROW>
end

if isempty(XYZ)
    error('G09_structure: no atoms read from step %d in %s.', step_idx, filename);
end

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
mol.symbols     = symbols;
mol.xyz         = XYZ;
mol.Z           = Zvec;
mol.Natoms      = size(XYZ, 1);
mol.step        = step_idx;
mol.n_steps     = n_blocks;
mol.orientation = 'Input orientation';
mol.filename    = filename;

end  % G09_structure


% =========================================================================
%  Local function: Z -> symbol table
% =========================================================================
function tbl = Z2symbol_table()
symbols = { ...
    'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', ...
    'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca', ...
    'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', ...
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr', ...
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', ...
    'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', ...
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', ...
    'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', ...
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', ...
    'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', ...
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', ...
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'};
keys_int = num2cell(int32(1:numel(symbols)));
tbl = containers.Map(keys_int, symbols);
end
