function hp = G16_hyperpolar(filename, varargin)
% G16_HYPERPOLAR  Extracts the dipole hyperpolarisability Beta from a file
%                 in Gaussian 16 .out/.log format.
%
%   hp = G16_HYPERPOLAR(filename)
%   hp = G16_HYPERPOLAR(filename, 'units', 'au')     % default
%   hp = G16_HYPERPOLAR(filename, 'units', 'esu')    % 10^-30 esu
%   hp = G16_HYPERPOLAR(filename, 'units', 'SI')     % 10^-50 C^3m^3J^-2
%
%   OUTPUT  struct hp with fields:
%   ── Static hyperpolarisability Beta(0;0,0) ────────────────────────────
%       .beta0              struct with the following fields:
%           .par_z          parallel component β‖ along z (au)
%           .perp_z         perpendicular component β⊥ w.r.t. z
%           .vec_x .vec_y .vec_z   vector components
%           .beta_vec       magnitude |β_vec| = sqrt(vec_x²+vec_y²+vec_z²)
%           .tensor         struct with all Cartesian components xxx,xxy,yxy,...,zzz
%           .units          stringa unita'
%
%   ── Dynamic hyperpolarisability Beta(-w;w,0) ──────────────────────────
%       .beta_dyn           struct array, one element per laser frequency:
%           .lambda_nm      laser wavelength (nm)
%           .par_z  .perp_z  .vec_x  .vec_y  .vec_z  .beta_vec
%           .tensor         struct with all Cartesian tensor components
%
%       .N_dyn              number of dynamic (laser) frequencies
%
%   ── Vibrational hyperpolarisability ────────────────────────────────
%       .beta_vib           [1x3] diagonal vibrational hyperpolarisability
%       .has_vib            logical        true if vibrational beta was found
%
%       .filename           char           source file path
%
%   Example:
%       hp = G16_hyperpolar('V_E00t.out');
%       hp.beta0.beta_vec          % modulo beta statico in au
%       hp.beta_dyn(1).lambda_nm   % wavelength of the first dynamic entry (nm)
%       hp.beta_dyn(1).beta_vec    % |beta| a quella frequenza

% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename', @ischar);
addParameter(p, 'units',    'au', @ischar);
parse(p, filename, varargin{:});
units = lower(p.Results.units);

if ~isfile(filename)
    error('G16_hyperpolar: file not found: %s', filename);
end
fid  = fopen(filename, 'r');
raw  = fread(fid, '*char')';
fclose(fid);
lines = strsplit(raw, newline);
N = numel(lines);

% -------------------------------------------------------------------------
% Conversions (G16 prints au, then 10^-30 esu, then 10^-50 SI)
% -------------------------------------------------------------------------
au2esu30 = 8.6392e-3;    % 1 au beta = 8.6392e-33 esu = 0.008639 x 10^-30 esu
au2SI50  = 3.2063e-3;    % 1 au beta = 3.2063e-53 SI = 0.003206 x 10^-50 SI

switch units
    case 'au',  ufac = 1;          ulbl = 'au';
    case 'esu', ufac = au2esu30;   ulbl = '10^{-30} esu';
    case 'si',  ufac = au2SI50;    ulbl = '10^{-50} C^3m^3J^{-2}';
    otherwise,  error('G16_hyperpolar: units deve essere ''au'', ''esu'' o ''SI''.');
end

% -------------------------------------------------------------------------
% Parsing
% -------------------------------------------------------------------------
beta0_raw   = [];    % raw (au) data for Beta(0;0,0)
beta_dyn_raw = struct('lambda_nm',{},'data',{});
beta_vib     = [];

k = 1;
while k <= N
    ln = lines{k};

    % ── Vibrational hyperpolarisability
    if ~isempty(strfind(ln, 'Diagonal vibrational hyperpolarizability:'))
        if k+1 <= N
            vals = sscanf(lines{k+1}, '%f');
            if numel(vals) >= 3
                beta_vib = vals(1:3)';
            end
        end
        k = k+2; continue
    end

    % ── Beta(0;0,0): static
    if ~isempty(regexp(ln, '^\s*Beta\(0;0,0\):', 'once'))
        beta0_raw = parse_beta_block(lines, k+1, N);
        k = k+15; continue
    end

    % ── Beta(-w;w,0): dynamic at frequency w
    tok_dyn = regexp(ln, 'Beta\(-w;w,0\)\s+w=\s*([\d.]+)\s*nm', 'tokens', 'once');
    if ~isempty(tok_dyn)
        lam = str2double(tok_dyn{1});
        bdata = parse_beta_block(lines, k+1, N);
        if ~isempty(bdata)
            nd = numel(beta_dyn_raw)+1;
            beta_dyn_raw(nd).lambda_nm = lam;
            beta_dyn_raw(nd).data      = bdata;
        end
        k = k+15; continue
    end

    k = k+1;
end

% -------------------------------------------------------------------------
% Assemble output struct and apply unit conversion
% -------------------------------------------------------------------------
hp.beta0   = apply_units(beta0_raw,   ufac, ulbl);
hp.beta_dyn = struct('lambda_nm',{},'par_z',{},'perp_z',{},...
                     'vec_x',{},'vec_y',{},'vec_z',{},...
                     'beta_vec',{},'tensor',{},'units',{});
for nd = 1:numel(beta_dyn_raw)
    entry = apply_units(beta_dyn_raw(nd).data, ufac, ulbl);
    entry.lambda_nm = beta_dyn_raw(nd).lambda_nm;
    hp.beta_dyn(nd) = entry;
end
hp.N_dyn = numel(beta_dyn_raw);

if ~isempty(beta_vib)
    hp.beta_vib  = beta_vib * ufac;
    hp.has_vib   = true;
else
    hp.beta_vib  = [];
    hp.has_vib   = false;
end
hp.filename = filename;

% -------------------------------------------------------------------------
% Print summary
% -------------------------------------------------------------------------
fprintf('\n── G16_hyperpolar: %s ──\n', filename);
if ~isempty(beta0_raw)
    fprintf('  Beta(0;0,0):  |β_vec| = %.2f %s\n', hp.beta0.beta_vec, ulbl);
    fprintf('    || (z) = %.2f   _|_(z) = %.2f\n', hp.beta0.par_z, hp.beta0.perp_z);
end
for nd = 1:hp.N_dyn
    fprintf('  Beta(-w;w,0) λ=%.1f nm:  |β_vec| = %.2f %s\n', ...
        hp.beta_dyn(nd).lambda_nm, hp.beta_dyn(nd).beta_vec, ulbl);
end
if hp.has_vib
    fprintf('  β_vib (diag) = [%.2f  %.2f  %.2f] %s\n', ...
        hp.beta_vib(1), hp.beta_vib(2), hp.beta_vib(3), ulbl);
end
fprintf('\n');

end  % G16_hyperpolar


% =========================================================================
%  Local function: parse Beta block
% =========================================================================
function out = parse_beta_block(lines, k_start, N)
out = [];

fields = struct('par_z',NaN,'perp_z',NaN,'vec_x',NaN,'vec_y',NaN,'vec_z',NaN,...
                'norm_par',NaN);
tens = struct('xxx',NaN,'xxy',NaN,'yxy',NaN,'yyy',NaN,...
              'xxz',NaN,'yxz',NaN,'yyz',NaN,'zxz',NaN,'zyz',NaN,'zzz',NaN,...
              'yxx',NaN,'zxx',NaN,'zyx',NaN,'zzx',NaN,...
              'xxy2',NaN,'yxy2',NaN,'zyy',NaN,'zzy',NaN,...
              'xxz2',NaN,'yxz2',NaN,'zzz2',NaN);

for k = k_start : min(k_start+25, N)
    ln = lines{k};
    % End of block: stop at blank lines or new section headers
    if ~isempty(regexp(ln, '^\s*$', 'once')), break; end
    if ~isempty(regexp(ln, '^\s*-{10,}', 'once')), break; end

    % Formato:  "   label   value_au   value_esu   value_SI"
    % The au value is the first number after the label
    tok = regexp(ln, '^\s*(\S+)\s+([-\d.Dd+E]+)', 'tokens', 'once');
    if isempty(tok), continue; end
    label = tok{1};
    val   = fortran2double(tok{2});

    switch label
        case '||(z)',   fields.par_z  = val;
        case '_|_(z)',  fields.perp_z = val;
        case '||',      fields.norm_par = val;
        case 'x',       fields.vec_x  = val;
        case 'y',       fields.vec_y  = val;
        case 'z',       fields.vec_z  = val;
        case 'xxx',  tens.xxx = val;
        case 'xxy',  tens.xxy = val;
        case 'yxy',  tens.yxy = val;
        case 'yyy',  tens.yyy = val;
        case 'xxz',  tens.xxz = val;
        case 'yxz',  tens.yxz = val;
        case 'yyz',  tens.yyz = val;
        case 'zxz',  tens.zxz = val;
        case 'zyz',  tens.zyz = val;
        case 'zzz',  tens.zzz = val;
        % Beta(-w;w,0) may list components in a different order
        case 'yxx',  tens.yxx = val;
        case 'zxx',  tens.zxx = val;
        case 'zyx',  tens.zyx = val;
        case 'zzx',  tens.zzx = val;
        case 'zyy',  tens.zyy = val;
        case 'zzy',  tens.zzy = val;
    end
end

if isnan(fields.vec_x), return; end

out.par_z    = fields.par_z;
out.perp_z   = fields.perp_z;
out.vec_x    = fields.vec_x;
out.vec_y    = fields.vec_y;
out.vec_z    = fields.vec_z;
out.beta_vec = sqrt(fields.vec_x^2 + fields.vec_y^2 + fields.vec_z^2);
out.tensor   = tens;
out.units    = 'au';
end


% =========================================================================
%  Local function: apply conversion factor
% =========================================================================
function out = apply_units(raw, ufac, ulbl)
if isempty(raw)
    out = struct('par_z',NaN,'perp_z',NaN,'vec_x',NaN,'vec_y',NaN,...
                 'vec_z',NaN,'beta_vec',NaN,'tensor',[],'units',ulbl);
    return
end
out = raw;
out.par_z    = raw.par_z    * ufac;
out.perp_z   = raw.perp_z   * ufac;
out.vec_x    = raw.vec_x    * ufac;
out.vec_y    = raw.vec_y    * ufac;
out.vec_z    = raw.vec_z    * ufac;
out.beta_vec = raw.beta_vec * ufac;
% Scale all tensor components
fnames = fieldnames(raw.tensor);
for fi = 1:numel(fnames)
    v = raw.tensor.(fnames{fi});
    if ~isnan(v)
        out.tensor.(fnames{fi}) = v * ufac;
    end
end
out.units = ulbl;
end


% =========================================================================
%  Local function: Fortran D/d exponent notation -> MATLAB double
% =========================================================================
function v = fortran2double(s)
v = str2double(strrep(strrep(s, 'D', 'e'), 'd', 'e'));
end
