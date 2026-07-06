function ch_out = G09_charges_fchk(mol, ch, varargin)
% G09_CHARGES_FCHK  Visualises atomic charges from a G09_fchk_read struct
%                   directly, without requiring a .log / .out file.
%
%   ch_out = G09_CHARGES_FCHK(mol, ch)
%   ch_out = G09_CHARGES_FCHK(mol, ch, Name, Value, ...)
%
%   Designed to work with the output of G09_fchk_read:
%
%       data   = G09_fchk_read('3typ.fchk');
%       ch_out = G09_charges_fchk(data.mol, data.ch);
%
%   The function mirrors the interface of G09_charges exactly, so the two
%   can be used interchangeably once data.mol and data.ch are available.
%
%   INPUTS:
%       mol   struct   geometry struct from G09_fchk_read (data.mol)
%                      Required fields: .symbols, .xyz, .Natoms, .filename
%
%       ch    struct   charge struct from G09_fchk_read (data.ch)
%                      Required fields: .charges, .symbols, .type, .Natoms
%                      Optional field:  .charges_H (H-summed, may be empty)
%
%   Optional parameters (Name-Value):
%       'mode'        - 'atom'  per-atom charges (default)
%                       'heavy' hydrogen charges summed onto heavy atoms
%                               (requires ch.charges_H to be non-empty;
%                               .fchk files do not contain H-summed data —
%                               use G09_charges on the .log file for that)
%       'plot'        - true (default) | false
%       'AtomScale'   - CPK sphere scale (default: 0.35)
%       'BondTol'     - bond detection tolerance (default: 1.30)
%       'FontSize'    - charge label font size in points (default: 8)
%       'ColorScale'  - 'RdBu'  red=positive / blue=negative (default)
%                       'none'  all labels in black
%       'threshold'   - hide labels where |q| < threshold (default: 0)
%
%   OUTPUT:
%       ch_out    struct with fields (same as G09_charges output):
%           .symbols    {Natoms x 1}   atomic symbols
%           .charges    [Natoms x 1]   per-atom charges (e)
%           .charges_H  [Nheavy x 1]  H-summed charges (empty if not available)
%           .sum_q      double         sum of all charges
%           .type       char           charge type, e.g. 'Mulliken'
%           .label      char           label string
%           .Natoms     int
%           .filename   char
%
%   Example:
%       data   = G09_fchk_read('3typ.fchk');
%
%       % Default: Mulliken per-atom charges with 3D visualisation
%       ch = G09_charges_fchk(data.mol, data.ch);
%
%       % Suppress plot, apply threshold
%       ch = G09_charges_fchk(data.mol, data.ch, 'plot', false, 'threshold', 0.05);
%
%       % Use as drop-in replacement for G09_charges:
%       ch = G09_charges_fchk(data.mol, data.ch, 'ColorScale', 'none');

% -------------------------------------------------------------------------
% Parse arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'mol',        @isstruct);
addRequired(p,  'ch',         @isstruct);
addParameter(p, 'mode',       'atom',   @ischar);
addParameter(p, 'plot',       true,     @islogical);
addParameter(p, 'AtomScale',  0.35,     @isnumeric);
addParameter(p, 'BondTol',    1.30,     @isnumeric);
addParameter(p, 'FontSize',   8,        @isnumeric);
addParameter(p, 'ColorScale', 'RdBu',   @ischar);
addParameter(p, 'threshold',  0,        @isnumeric);
parse(p, mol, ch, varargin{:});

mode       = lower(p.Results.mode);
do_plot    = p.Results.plot;
atom_scale = p.Results.AtomScale;
bond_tol   = p.Results.BondTol;
fsize      = p.Results.FontSize;
cscale     = p.Results.ColorScale;
thr        = p.Results.threshold;

% -------------------------------------------------------------------------
% Validate inputs
% -------------------------------------------------------------------------
required_mol = {'symbols', 'xyz', 'Natoms'};
for k = 1:numel(required_mol)
    if ~isfield(mol, required_mol{k})
        error('G09_charges_fchk: mol is missing field "%s". Use data.mol from G09_fchk_read.', ...
            required_mol{k});
    end
end

required_ch = {'charges', 'symbols', 'type', 'Natoms'};
for k = 1:numel(required_ch)
    if ~isfield(ch, required_ch{k})
        error('G09_charges_fchk: ch is missing field "%s". Use data.ch from G09_fchk_read.', ...
            required_ch{k});
    end
end

if mol.Natoms ~= ch.Natoms
    error('G09_charges_fchk: mol.Natoms (%d) ≠ ch.Natoms (%d).', ...
        mol.Natoms, ch.Natoms);
end

% -------------------------------------------------------------------------
% Select charge set according to mode
% -------------------------------------------------------------------------
q_atom  = ch.charges(:);
q_heavy = [];
if isfield(ch, 'charges_H') && ~isempty(ch.charges_H)
    q_heavy = ch.charges_H(:);
end

switch mode
    case 'atom'
        syms_use = mol.symbols;
        xyz_use  = mol.xyz;
        q_use    = q_atom;

    case 'heavy'
        if isempty(q_heavy)
            warning(['G09_charges_fchk: H-summed charges (charges_H) are not available ' ...
                'in .fchk files. Falling back to per-atom charges.\n' ...
                'For H-summed charges, use G09_charges on the .log file.']);
            syms_use = mol.symbols;
            xyz_use  = mol.xyz;
            q_use    = q_atom;
        else
            is_heavy = ~strcmp(mol.symbols, 'H');
            syms_use = mol.symbols(is_heavy);
            xyz_use  = mol.xyz(is_heavy, :);
            q_use    = q_heavy;
            if numel(q_use) ~= size(xyz_use, 1)
                % Fallback: take first numel(q_use) heavy atoms
                xyz_use = xyz_use(1:numel(q_use), :);
            end
        end

    otherwise
        error('G09_charges_fchk: mode must be ''atom'' or ''heavy''.');
end

% -------------------------------------------------------------------------
% Build output struct  (identical layout to G09_charges)
% -------------------------------------------------------------------------
src = '';
if isfield(mol, 'filename'), src = mol.filename; end
if isfield(ch, 'filename'),  src = ch.filename;  end

charge_type  = ch.type;
found_label  = '';
if isfield(ch, 'label'), found_label = ch.label; end
if isempty(found_label)
    found_label = sprintf('%s Charges (from .fchk)', charge_type);
end

ch_out.symbols   = mol.symbols(:);
ch_out.charges   = q_atom;
ch_out.charges_H = q_heavy;
ch_out.sum_q     = sum(q_atom);
ch_out.type      = charge_type;
ch_out.label     = found_label;
ch_out.Natoms    = mol.Natoms;
ch_out.filename  = src;

% -------------------------------------------------------------------------
% Print table
% -------------------------------------------------------------------------
[~, fname] = fileparts(src);
fprintf('\n── G09_charges_fchk (%s, %s): %s ──\n', charge_type, mode, fname);
fprintf('  Source : %s\n', found_label);
fprintf('  %4s  %-4s  %8s\n', 'Idx', 'Sym', 'q (e)');
fprintf('  %s\n', repmat('-', 1, 22));
for i = 1:numel(q_use)
    fprintf('  %4d  %-4s  %+8.4f\n', i, syms_use{i}, q_use(i));
end
fprintf('  %s\n', repmat('-', 1, 22));
fprintf('  Sum = %+.5f e\n\n', sum(q_use));

% -------------------------------------------------------------------------
% 3D visualisation
% -------------------------------------------------------------------------
if ~do_plot
    return
end

fig = figure('Color', 'white', ...
    'Name',  sprintf('%s charges — %s', charge_type, fname), ...
    'NumberTitle', 'off');
ax = axes('Parent', fig);

G09_draw_molecule(mol, 'Ax', ax, ...
    'AtomScale',  atom_scale, ...
    'BondTol',    bond_tol, ...
    'ShowLabels', false, ...
    'ShowLegend', true, ...
    'Title', sprintf('%s — %s charges (%s)', ...
        strrep(fname, '_', '\_'), charge_type, mode));

q_max = max(abs(q_use));
if q_max == 0, q_max = 1; end

hold(ax, 'on');
for i = 1:numel(q_use)
    if abs(q_use(i)) < thr, continue; end

    if strcmpi(cscale, 'RdBu')
        clr = charge_color(q_use(i) / q_max);
    else
        clr = [0.05 0.05 0.05];
    end

    r_off = atom_scale * 0.8 + 0.3;
    text(ax, xyz_use(i,1), xyz_use(i,2), xyz_use(i,3) + r_off, ...
         sprintf('%+.3f', q_use(i)), ...
         'FontSize',            fsize, ...
         'Color',               clr, ...
         'FontWeight',          'bold', ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment',   'bottom', ...
         'Interpreter',         'none', ...
         'HandleVisibility',    'off');
end

rotate3d(ax, 'on');

end  % G09_charges_fchk


% =========================================================================
%  RdBu colour map: red (positive) – white (zero) – blue (negative)
% =========================================================================
function clr = charge_color(t)
% t in [-1, +1]
t = max(-1, min(1, t));
if t >= 0                        % positive → red
    clr = [1, 1-t, 1-t];
else                             % negative → blue
    clr = [1+t, 1+t, 1];
end
clr = clr * 0.82;
end