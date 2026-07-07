function ch = G16_charges(filename, varargin)
% G16_CHARGES  Extracts Mulliken or APT atomic charges from a Gaussian 16
%              output file and optionally renders them on the 3D structure.
%
%   ch = G16_CHARGES(filename)
%   ch = G16_CHARGES(filename, 'type', 'APT')
%   ch = G16_CHARGES(filename, 'mode', 'heavy')
%   ch = G16_CHARGES(filename, 'plot', true, 'threshold', 0.05)
%
%   ROBUSTNESS: Uses fuzzy label matching to handle all known charge header
%   variants across Gaussian 16 revisions:
%
%     G16 (all revisions):       "Mulliken charges:"
%     G09-compatible builds:     "Mulliken atomic charges:"
%
%   The parser finds any line containing the charge type keyword and
%   "charges:", excluding "Sum of" and H-summed lines automatically.
%   The exact label found in the file is reported in ch.label.
%
%   Optional parameters (Name-Value):
%       'type'            - 'Mulliken' (default) | 'APT'
%       'mode'            - 'atom' (default) | 'heavy'
%       'plot'            - true  renders labels on 3D structure (default: true)
%       'AtomScale'       - CPK sphere scale (default: 0.35)
%       'BondTol'         - bond detection tolerance (default: 1.30)
%       'FontSize'        - label font size in points (default: 8)
%       'ColorScale'      - 'RdBu' blue=neg/red=pos (default) | 'none'
%       'threshold'       - hide labels with |q| < threshold (default: 0)
%       'ShowDipole'      - overlay the total dipole moment vector,
%                           anchored at the charge-weighted centroid of
%                           the negative ("electronic") partial charges
%                           by default (default: false)
%       'DipoleOrigin'    - anchor point for the dipole arrow:
%                             'negcharge' (default) - centroid of
%                                 atoms with NEGATIVE partial charge
%                                 (the "center of electronic charge")
%                             'poscharge'  - centroid of positive charges
%                             'centroid'   - unweighted atom centroid
%                             [x y z]      - explicit point
%       'DipoleScale'     - arrow length, Angstrom per Debye (default: 1.0)
%       'DipoleColor'     - arrow colour (default: [0 0.6 0])
%       'DipoleLineWidth' - arrow line width (default: 2.5)
%       'ShowDipoleLabel' - annotate the arrow with |mu| in Debye (default: true)
%       'DipoleFontSize'  - dipole label font size (default: 11)
%       'DipoleUnits'     - units used to display |mu| (arrow label and
%                           command-window table) when 'ShowDipole' is
%                           true: 'Debye' (default) | 'au'
%
%   OUTPUT  struct ch with fields:
%       .symbols        {Natoms x 1}   atomic symbols
%       .charges        [Natoms x 1]   per-atom charges (e)
%       .charges_H      [Nheavy x 1]   H-summed charges on heavy atoms
%       .sum_q          double         sum of all charges (≈ 0 for neutral)
%       .type           char           'Mulliken' or 'APT'
%       .label          char           exact header label found in file
%       .Natoms         int
%       .filename       char
%       .dipole         [1x3]          dipole moment (Debye), only when
%                                      'ShowDipole' is true; [] otherwise
%       .dipole_origin  [1x3]          anchor point used for the arrow
%                                      (Angstrom), only when 'ShowDipole'
%                                      is true; [] otherwise
%       .dipole_Debye   double         |mu| in Debye, only when
%                                      'ShowDipole' is true; [] otherwise
%       .dipole_au      double         |mu| in atomic units, only when
%                                      'ShowDipole' is true; [] otherwise

% -------------------------------------------------------------------------
% Parse arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename',   @ischar);
addParameter(p, 'type',       'Mulliken', @ischar);
addParameter(p, 'mode',       'atom',     @ischar);
addParameter(p, 'plot',       true,       @islogical);
addParameter(p, 'AtomScale',  0.35,       @isnumeric);
addParameter(p, 'BondTol',    1.30,       @isnumeric);
addParameter(p, 'FontSize',   8,          @isnumeric);
addParameter(p, 'ColorScale', 'RdBu',     @ischar);
addParameter(p, 'threshold',  0,          @isnumeric);
addParameter(p, 'ShowDipole',      false,      @islogical);
addParameter(p, 'DipoleOrigin',    'negcharge',@(x) ischar(x) || (isnumeric(x) && numel(x)==3));
addParameter(p, 'DipoleScale',     1.0,        @isnumeric);
addParameter(p, 'DipoleColor',     [0 0.6 0],  @(x) isnumeric(x) && numel(x)==3);
addParameter(p, 'DipoleLineWidth', 2.5,        @isnumeric);
addParameter(p, 'ShowDipoleLabel', true,       @islogical);
addParameter(p, 'DipoleFontSize',  11,         @isnumeric);
addParameter(p, 'DipoleUnits',     'Debye',    @ischar);
addParameter(p, 'Lines',           {},         @iscell);
parse(p, filename, varargin{:});

charge_type = p.Results.type;
mode        = lower(p.Results.mode);
do_plot     = p.Results.plot;
atom_scale  = p.Results.AtomScale;
bond_tol    = p.Results.BondTol;
fsize       = p.Results.FontSize;
cscale      = p.Results.ColorScale;
thr         = p.Results.threshold;
show_dipole  = p.Results.ShowDipole;
dip_origin   = p.Results.DipoleOrigin;
dip_scale    = p.Results.DipoleScale;
dip_color    = p.Results.DipoleColor;
dip_lw       = p.Results.DipoleLineWidth;
dip_label_on = p.Results.ShowDipoleLabel;
dip_fsize    = p.Results.DipoleFontSize;
dip_units    = p.Results.DipoleUnits;

switch lower(dip_units)
    case 'debye'
        dip_units = 'Debye';
    case 'au'
        dip_units = 'au';
    otherwise
        warning('G16_charges:badDipoleUnits', ...
            'Unknown ''DipoleUnits'' = "%s"; using ''Debye''.', dip_units);
        dip_units = 'Debye';
end
DEBYE_TO_AU = 0.393430;   % 1 Debye = 0.393430 a.u. (1 a.u. = 2.541746 Debye)

% -------------------------------------------------------------------------
% Read file
% -------------------------------------------------------------------------
lines = p.Results.Lines;
if isempty(lines)
    if ~isfile(filename)
        error('G16_charges: file not found: %s', filename);
    end
    fid  = fopen(filename, 'r');
    raw  = fread(fid, '*char')';
    fclose(fid);
    lines = strsplit(raw, newline);
end
N = numel(lines);

% -------------------------------------------------------------------------
% Locate charge blocks using fuzzy matching
% -------------------------------------------------------------------------
% Finds any line that:
%   (a) contains the charge type keyword (case-insensitive)
%   (b) contains "charges"
%   (c) ends with ':'
%   (d) does NOT contain "Sum of"
% Then classifies as per-atom vs H-summed based on "hydrogen"/"summed".

atom_starts  = [];
heavy_starts = [];
found_labels = {};

for k = 1:N
    ln     = lines{k};
    ln_low = lower(ln);

    if ~contains(ln_low, lower(charge_type)), continue; end
    if ~contains(ln_low, 'charges'),          continue; end

    ln_trim = strtrim(ln);
    if isempty(ln_trim) || ln_trim(end) ~= ':', continue; end

    if contains(ln_low, 'sum of'), continue; end

    if contains(ln_low, 'hydrogen') || contains(ln_low, 'summed')
        heavy_starts(end+1) = k; %#ok<AGROW>
    else
        atom_starts(end+1) = k;          %#ok<AGROW>
        found_labels{end+1} = ln_trim;   %#ok<AGROW>
    end
end

% -------------------------------------------------------------------------
% Validate
% -------------------------------------------------------------------------
if isempty(atom_starts)
    candidates = {};
    for k = 1:min(N, 300)
        if contains(lower(lines{k}), lower(charge_type))
            candidates{end+1} = sprintf('  line %d: %s', k, strtrim(lines{k})); %#ok<AGROW>
        end
    end
    if isempty(candidates)
        error('G16_charges: no "%s" charges found in %s', charge_type, filename);
    else
        error('G16_charges: charge header not found in %s\n"%s" appears in:\n%s', ...
            filename, charge_type, strjoin(candidates, '\n'));
    end
end

% -------------------------------------------------------------------------
% Helper: parse a charge block (per-atom or H-summed)
% -------------------------------------------------------------------------
    function [syms_out, q_out] = parse_block(k_start)
        syms_out = {};
        q_out    = [];
        k2 = k_start + 2;      % skip index row "               1"
        while k2 <= N
            ln2 = strtrim(lines{k2});
            if isempty(ln2),                              break; end
            if contains(lower(ln2), 'sum of'),            break; end
            if contains(lower(ln2), 'charges'),           break; end
            m = regexp(lines{k2}, ...
                '^\s*\d+\s+([A-Za-z]+)\s+([-\d.]+)', 'tokens', 'once');
            if ~isempty(m)
                syms_out{end+1} = m{1};             %#ok<AGROW>
                q_out(end+1)    = str2double(m{2}); %#ok<AGROW>
            end
            k2 = k2 + 1;
        end
        q_out = q_out(:);
    end

% -------------------------------------------------------------------------
% Parse per-atom block (last occurrence)
% -------------------------------------------------------------------------
[syms_atom, q_atom] = parse_block(atom_starts(end));

if isempty(q_atom)
    error('G16_charges: charge header found but no atom data read from %s', filename);
end

found_label = found_labels{end};

% -------------------------------------------------------------------------
% Parse H-summed block (last occurrence, if present)
% -------------------------------------------------------------------------
syms_heavy = {};
q_heavy    = [];
if ~isempty(heavy_starts)
    [syms_heavy, q_heavy] = parse_block(heavy_starts(end));
end

% -------------------------------------------------------------------------
% Select data set according to mode
% -------------------------------------------------------------------------
switch mode
    case 'atom'
        syms_use = syms_atom;
        q_use    = q_atom;
    case 'heavy'
        if isempty(q_heavy)
            warning('G16_charges: H-summed charges not found; falling back to per-atom.');
            syms_use = syms_atom;
            q_use    = q_atom;
        else
            syms_use = syms_heavy;
            q_use    = q_heavy;
        end
    otherwise
        error('G16_charges: mode must be ''atom'' or ''heavy''.');
end

% -------------------------------------------------------------------------
% Build output struct
% -------------------------------------------------------------------------
ch.symbols   = syms_atom(:);
ch.charges   = q_atom;
ch.charges_H = q_heavy;
ch.sum_q     = sum(q_atom);
ch.type      = charge_type;
ch.label     = found_label;
ch.Natoms    = numel(q_atom);
ch.filename  = filename;
ch.dipole        = [];
ch.dipole_origin = [];
ch.dipole_Debye  = [];
ch.dipole_au     = [];

% -------------------------------------------------------------------------
% Dipole moment (optional). Computed here -- rather than inside the
% plotting block below -- so that .dipole / .dipole_origin are populated
% in the output struct even when 'plot' is false.
% -------------------------------------------------------------------------
mol = [];   % loaded lazily, at most once, by whichever block needs it first
if show_dipole
    try
        dp = G16_dipole_polar(filename, 'units', 'Debye', 'Lines', lines);
        mu = local_extract_dipole(dp);
    catch ME
        warning('G16_charges:dipoleReadFailed', ...
            'Could not read the dipole moment (%s); ''ShowDipole'' will be ignored.', ME.message);
        mu = [];
    end
    if isempty(mu) || norm(mu) < eps
        if isempty(mu)
            warning('G16_charges:dipoleFieldNotFound', ...
                'Dipole moment field not recognised in G16_dipole_polar output; ''ShowDipole'' will be ignored.');
        end
    else
        mol = G16_structure(filename);
        ch.dipole        = mu;
        ch.dipole_origin = local_dipole_origin(dip_origin, mol.xyz, q_atom);
        ch.dipole_Debye  = norm(mu);
        ch.dipole_au     = norm(mu) * DEBYE_TO_AU;
    end
end

% Value and unit symbol used for the dipole magnitude wherever it is
% displayed (command-window table and 3D arrow label), per 'DipoleUnits'.
if ~isempty(ch.dipole)
    if strcmp(dip_units, 'au')
        dip_value_disp  = ch.dipole_au;
        dip_unit_symbol = 'a.u.';
    else
        dip_value_disp  = ch.dipole_Debye;
        dip_unit_symbol = 'D';
    end
end

% Print table
fprintf('\n── G16_charges (%s, %s): %s ──\n', charge_type, mode, filename);
fprintf('  Header found : "%s"\n', found_label);
fprintf('  %4s  %-4s  %8s\n', 'Idx', 'Sym', 'q (e)');
fprintf('  %s\n', repmat('-', 1, 22));
for i = 1:numel(syms_use)
    fprintf('  %4d  %-4s  %+8.4f\n', i, syms_use{i}, q_use(i));
end
fprintf('  %s\n', repmat('-', 1, 22));
fprintf('  Sum = %+.5f e\n', sum(q_use));
if ~isempty(ch.dipole)
    fprintf('  |mu| = %.3f %s  (anchor: %s)\n', dip_value_disp, dip_unit_symbol, local_origin_label(dip_origin));
end
fprintf('\n');

% -------------------------------------------------------------------------
% 3D visualisation
% -------------------------------------------------------------------------
if do_plot
    if isempty(mol)
        mol = G16_structure(filename);
    end
    [~, fname] = fileparts(filename);
    fig = figure('Color', 'white', ...
        'Name', sprintf('%s charges', charge_type), 'NumberTitle', 'off');
    ax  = axes('Parent', fig);
    G16_draw_molecule(mol, 'Ax', ax, 'AtomScale', atom_scale, ...
        'BondTol', bond_tol, 'ShowLabels', false, 'ShowLegend', true, ...
        'Title', sprintf('%s — %s charges (%s)', ...
            strrep(fname, '_', '\_'), charge_type, mode));

    % Coordinates for label placement
    if strcmp(mode, 'atom')
        xyz_use = mol.xyz;
    else
        is_heavy = ~strcmp(mol.symbols, 'H');
        xyz_use  = mol.xyz(is_heavy, :);
        if size(xyz_use, 1) ~= numel(q_use)
            xyz_use = mol.xyz(1:numel(q_use), :);
        end
    end

    q_max = max(abs(q_use));
    if q_max == 0, q_max = 1; end

    hold(ax, 'on');
    for i = 1:numel(q_use)
        if abs(q_use(i)) < thr, continue; end
        if strcmpi(cscale, 'RdBu')
            clr = charge_color(q_use(i) / q_max);
        else
            clr = [0 0 0];
        end
        r_off = atom_scale * 0.8 + 0.3;
        text(ax, xyz_use(i,1), xyz_use(i,2), xyz_use(i,3) + r_off, ...
             sprintf('%+.3f', q_use(i)), ...
             'FontSize', fsize, 'Color', clr, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'Interpreter', 'none', 'HandleVisibility', 'off');
    end
    rotate3d(ax, 'on');

    % ---------------------------------------------------------------
    % Dipole moment arrow (optional)
    % ---------------------------------------------------------------
    if show_dipole && ~isempty(ch.dipole)
        mu     = ch.dipole;
        origin = ch.dipole_origin;
        vec    = (mu / norm(mu)) * norm(mu) * dip_scale;   % cosmetic length only

        quiver3(ax, origin(1), origin(2), origin(3), vec(1), vec(2), vec(3), 0, ...
            'Color', dip_color, 'LineWidth', dip_lw, 'MaxHeadSize', 0.6, ...
            'HandleVisibility', 'off');

        if dip_label_on
            tip = origin + vec;
            text(ax, tip(1), tip(2), tip(3), sprintf('  \\mu = %.2f %s', dip_value_disp, dip_unit_symbol), ...
                 'Color', dip_color, 'FontSize', dip_fsize, 'FontWeight', 'bold', ...
                 'Interpreter', 'none', 'HandleVisibility', 'off');
        end
    end
end

end  % G16_charges


% =========================================================================
%  RdBu colour map: blue (negative) -> white (zero) -> red (positive)
% =========================================================================
function clr = charge_color(t)
t = max(-1, min(1, t));
if t >= 0
    clr = [1, 1-t, 1-t];
else
    clr = [1+t, 1+t, 1];
end
clr = clr * 0.82;
end


% =========================================================================
%  Dipole moment helpers
% =========================================================================
function mu = local_extract_dipole(dp)
%LOCAL_EXTRACT_DIPOLE  Pulls a [1x3] dipole vector (Debye) out of whatever
%   struct/array G16_dipole_polar returns, without assuming one exact
%   field name.
mu = [];
if isnumeric(dp) && numel(dp) == 3
    mu = double(dp(:))';
    return
end
if ~isstruct(dp), return; end

candidates = {'dipole', 'Dipole', 'mu', 'Mu', 'dip', 'vector'};
for i = 1:numel(candidates)
    if isfield(dp, candidates{i})
        v = dp.(candidates{i});
        if isnumeric(v) && numel(v) == 3
            mu = double(v(:))';
            return
        end
    end
end
if isfield(dp, 'mu_x') && isfield(dp, 'mu_y') && isfield(dp, 'mu_z')
    mu = double([dp.mu_x, dp.mu_y, dp.mu_z]);
end
end

% -------------------------------------------------------------------------
function origin = local_dipole_origin(originSpec, xyz, charges)
%LOCAL_DIPOLE_ORIGIN  Resolves the 'DipoleOrigin' option to a 3D point.
if isnumeric(originSpec) && numel(originSpec) == 3
    origin = originSpec(:)';
    return
end
switch lower(originSpec)
    case 'negcharge'
        origin = local_charge_centroid(xyz, charges, charges < 0, 'negative');
    case 'poscharge'
        origin = local_charge_centroid(xyz, charges, charges > 0, 'positive');
    case 'centroid'
        origin = mean(xyz, 1);
    otherwise
        warning('G16_charges:badDipoleOrigin', ...
            'Unknown ''DipoleOrigin'' = "%s"; using the unweighted atom centroid.', originSpec);
        origin = mean(xyz, 1);
end
end

% -------------------------------------------------------------------------
function origin = local_charge_centroid(xyz, charges, mask, label)
if ~any(mask)
    warning('G16_charges:noMatchingCharges', ...
        'No atoms with a %s partial charge; using the unweighted atom centroid.', label);
    origin = mean(xyz, 1);
    return
end
w = abs(charges(mask));
origin = sum(xyz(mask, :) .* w, 1) / sum(w);
end

% -------------------------------------------------------------------------
function s = local_origin_label(originSpec)
if isnumeric(originSpec)
    s = mat2str(originSpec, 3);
else
    s = originSpec;
end
end