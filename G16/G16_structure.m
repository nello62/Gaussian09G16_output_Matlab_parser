function mol = G16_structure(filename, varargin)
% G16_STRUCTURE  Extracts the molecular geometry from a Gaussian 16 .out/.log file.
%
%   mol = G16_STRUCTURE(filename)
%   mol = G16_STRUCTURE(filename, 'orientation', 'standard')
%   mol = G16_STRUCTURE(filename, 'orientation', 'input')
%   mol = G16_STRUCTURE(filename, 'step', N)
%   mol = G16_STRUCTURE(filename, 'step', 'last')   % default
%   mol = G16_STRUCTURE(filename, 'step', 'first')
%
% OUTPUT  struct mol with fields:
%   .symbols   {Natoms×1 cell}   atomic symbols ('C','H','N',...)
%   .xyz       [Natoms×3 double] Cartesian coordinates in Angstrom  (X Y Z)
%   .Z         [Natoms×1 int]    atomic numbers
%   .Natoms    int               number of atoms
%   .step      int               index of the extracted step (1-based)
%   .orientation  char           'standard' or 'input'
%   .filename  char              source file path
%
% EXAMPLES:
%   mol = G16_structure('zeatin.out');
%   mol = G16_structure('opt_traj.out', 'step', 1);   % initial geometry
%   mol = G16_structure('calc.out', 'orientation', 'input');
%
% NOTES:
%   - For opt+freq concatenated jobs, the last geometry block is the optimised one.
%     the geometry is the optimised one.
%   - If Standard orientation is absent (e.g. NoSymm keyword), Input orientation is used.
%   - step='last' (default) always takes the last geometry block in the file.

% -------------------------------------------------------------------------
% Parse input arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename',    @ischar);
addParameter(p, 'orientation', 'auto',  @(x) ischar(x) && ...
    any(strcmpi(x, {'standard','input','auto'})));
addParameter(p, 'step',        'last',  @(x) ischar(x) || isnumeric(x));
addParameter(p, 'Lines',       {},      @iscell);
parse(p, filename, varargin{:});

ori_pref = lower(p.Results.orientation);
step_req = p.Results.step;

% -------------------------------------------------------------------------
% Read file
% -------------------------------------------------------------------------
lines = p.Results.Lines;
if isempty(lines)
    if ~isfile(filename)
        error('G16_structure: file not found: %s', filename);
    end
    fid = fopen(filename, 'r');
    raw = fread(fid, '*char')';
    fclose(fid);
    lines = strsplit(raw, newline);
    G16_check_gaussian_match(lines, filename);
end

% -------------------------------------------------------------------------
% Atomic number -> symbol table
% -------------------------------------------------------------------------
Z2sym = Z2symbol_table();

% -------------------------------------------------------------------------
% Find all orientation blocks in the file
% -------------------------------------------------------------------------
% Search for "Standard orientation:" or "Input orientation:"
% Each block has the following structure:
%
%  Standard orientation:
%  ----------------------------------------------------
%   Center     Atomic      Atomic             Coordinates (Angstroms)
%   Number     Number       Type             X           Y           Z
%  ----------------------------------------------------
%       1          6           0        0.000000    0.000000    0.123456
%       ...
%  ----------------------------------------------------
%
% The block ends with the dashes line after the data.

% Orientation priority
if strcmp(ori_pref, 'auto')
    % use Standard if present, otherwise Input
    has_std = any(~cellfun(@isempty, regexpi(lines, 'Standard orientation\s*:')));
    if has_std
        ori_label = 'Standard orientation';
    else
        ori_label = 'Input orientation';
    end
elseif strcmp(ori_pref, 'standard')
    ori_label = 'Standard orientation';
else
    ori_label = 'Input orientation';
end

% Find the header lines of all steps
block_starts = find(~cellfun(@isempty, regexpi(lines, [ori_label, '\s*:'])));

if isempty(block_starts)
    % Fallback: if Standard was requested but not found, try Input
    if strcmp(ori_pref, 'standard')
        warning('G16_structure: Standard orientation not found, falling back to Input orientation.');
        ori_label   = 'Input orientation';
        block_starts = find(~cellfun(@isempty, regexpi(lines, [ori_label, '\s*:'])));
    end
    if isempty(block_starts)
        error('G16_structure: no "%s" block found in %s', ori_label, filename);
    end
end

% Step selection
n_blocks = numel(block_starts);
if ischar(step_req)
    switch lower(step_req)
        case 'last',  step_idx = n_blocks;
        case 'first', step_idx = 1;
        otherwise,    error('G16_structure: step must be ''first'', ''last'', or an integer.');
    end
else
    step_idx = round(step_req);
    if step_idx < 1 || step_idx > n_blocks
        error('G16_structure: step %d out of range [1, %d].', step_idx, n_blocks);
    end
end

header_line = block_starts(step_idx);

% -------------------------------------------------------------------------
% Parse the selected block
% -------------------------------------------------------------------------
% Skip the 5 header lines:
%   +1  separator  (-----)
%   +2  "Center  Atomic  Atomic  Coordinates (Angstroms)"
%   +3  "Number  Number  Type    X  Y  Z"
%   +4  separator  (-----)
%   +5  first data line
data_start = header_line + 5;

symbols = {};
XYZ     = [];
Zvec    = [];

for k = data_start : numel(lines)
    ln = strtrim(lines{k});
    if isempty(ln), continue; end
    
    % End of block: dashes-only line (after strtrim)
    if ~isempty(ln) && all(ln == '-'), break; end
    
    % Parse columns: CenterNum  AtomicNum  AtomType  X  Y  Z
    % Typical format:  "    1          6           0        0.000000  ..."
    tok = sscanf(ln, '%d %d %d %f %f %f');
    if numel(tok) ~= 6, continue; end   % skip non-data lines (e.g. empty lines)
    
    znum = tok(2);
    x    = tok(4);
    y    = tok(5);
    z    = tok(6);
    
    sym = Z2sym(znum);
    symbols{end+1, 1} = sym;           %#ok<AGROW>
    XYZ(end+1, :)     = [x, y, z];    %#ok<AGROW>
    Zvec(end+1, 1)    = znum;          %#ok<AGROW>
end

if isempty(XYZ)
    error('G16_structure: no atoms read from step %d.', step_idx);
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
mol.orientation = ori_label;
mol.filename    = filename;

end  % G16_structure


% =========================================================================
%  Local function: Z -> symbol table (Z = 1..118)
% =========================================================================
function tbl = Z2symbol_table()
% Returns a containers.Map: atomic number (int) -> symbol (char)

symbols = { ...
    'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', ...   1-10
    'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca', ...  11-20
    'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', ...  21-30
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr', ...  31-40
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', ...  41-50
    'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', ...  51-60
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', ...  61-70
    'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', ...  71-80
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', ...  81-90
    'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', ...  91-100
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', ... 101-110
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'};               % 111-118

keys_int   = num2cell(int32(1:numel(symbols)));
tbl = containers.Map(keys_int, symbols);
end
