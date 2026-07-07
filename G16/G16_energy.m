function en = G16_energy(filename, varargin)
% G16_ENERGY  Extracts energies from a Gaussian 16 .out/.log file.
%
%   en = G16_ENERGY(filename)
%   en = G16_ENERGY(filename, 'step', 'last')   % default
%   en = G16_ENERGY(filename, 'step', 'first')
%   en = G16_ENERGY(filename, 'step', N)        % N-esima occorrenza SCF
%
%   OUTPUT  struct en with fields (all values in Hartree unless stated otherwise):
%       .SCF          SCF energy (Hartree)
%       .method       method string  es. 'RB3LYP'
%       .ZPE_corr     zero-point energy correction (Hartree)
%       .U_corr       thermal correction to energy (Hartree)
%       .H_corr       thermal correction to enthalpy (Hartree)
%       .G_corr       thermal correction to Gibbs free energy (Hartree)
%       .E0           SCF + ZPE  (Hartree)
%       .U            SCF + U_corr (Hartree)
%       .H            SCF + H_corr (Hartree)
%       .G            SCF + G_corr (Hartree)
%       .ZPE_kJ       Zero-point vibrational energy (kJ/mol)
%       .T            Temperature (K)
%       .P            Pressure (atm)
%       .has_thermo   logical — true if the thermochemistry section was found
%       .filename     char
%
%   Note: for opt-only jobs (no freq), has_thermo = false
%         e the *_corr, E0, U, H, G, T, P fields will be NaN.

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
% Read file
% -------------------------------------------------------------------------
lines = p.Results.Lines;
if isempty(lines)
    if ~isfile(filename)
        error('G16_energy: file not found: %s', filename);
    end
    fid  = fopen(filename,'r');
    raw  = fread(fid,'*char')';
    fclose(fid);
    lines = strsplit(raw, newline);
end

% -------------------------------------------------------------------------
% SCF Done  (may appear multiple times)
% -------------------------------------------------------------------------
scf_vals   = [];
scf_methods = {};

for k = 1:numel(lines)
    ln = lines{k};
    % SCF Done:  E(RB3LYP) =  -1160.20272617     A.U. after ...
    tok = regexp(ln, 'SCF Done:\s+E\((\S+)\)\s*=\s*([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok)
        scf_methods{end+1} = tok{1};   %#ok<AGROW>
        scf_vals(end+1)    = str2double(tok{2}); %#ok<AGROW>
    end
end

if isempty(scf_vals)
    error('G16_energy: no "SCF Done" line found in %s', filename);
end

% Step selection
n_scf = numel(scf_vals);
if ischar(step_req)
    switch lower(step_req)
        case 'last',  si = n_scf;
        case 'first', si = 1;
        otherwise, error('G16_energy: step must be ''first'', ''last'' or an integer.');
    end
else
    si = round(step_req);
    if si < 1 || si > n_scf
        error('G16_energy: step %d out of range [1,%d].', si, n_scf);
    end
end

en.SCF    = scf_vals(si);
en.method = scf_methods{si};

% -------------------------------------------------------------------------
% Thermal corrections  (takes the LAST occurrence in file — freq section)
% -------------------------------------------------------------------------
% Patterns for thermochemistry lines:
%  Zero-point correction=                           0.291821 (Hartree/Particle)
%  Thermal correction to Energy=                    0.311936
%  Thermal correction to Enthalpy=                  0.312880
%  Thermal correction to Gibbs Free Energy=         0.242309
%  Sum of electronic and zero-point Energies=      -1159.910939
%  Sum of electronic and thermal Energies=         -1159.890825
%  Sum of electronic and thermal Enthalpies=       -1159.889880
%  Sum of electronic and thermal Free Energies=    -1159.960451

ZPE_corr = NaN; U_corr = NaN; H_corr = NaN; G_corr = NaN;
E0 = NaN; U_tot = NaN; H_tot = NaN; G_tot = NaN;
ZPE_kJ = NaN; T = NaN; P_atm = NaN;

for k = 1:numel(lines)
    ln = lines{k};

    % Temperature e Pressure
    tok = regexp(ln, 'Temperature\s+([\d.]+)\s+Kelvin\.\s+Pressure\s+([\d.]+)', 'tokens', 'once');
    if ~isempty(tok)
        T   = str2double(tok{1});
        P_atm = str2double(tok{2});
        continue
    end

    % ZPE in J/mol
    tok = regexp(ln, 'Zero-point vibrational energy\s+([\d.]+)\s+\(Joules', 'tokens', 'once');
    if ~isempty(tok)
        ZPE_kJ = str2double(tok{1}) / 1000;   % J/mol -> kJ/mol
        continue
    end

    % Correzioni (in Hartree)
    tok = regexp(ln, 'Zero-point correction=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), ZPE_corr = str2double(tok{1}); continue; end

    tok = regexp(ln, 'Thermal correction to Energy=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), U_corr = str2double(tok{1}); continue; end

    tok = regexp(ln, 'Thermal correction to Enthalpy=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), H_corr = str2double(tok{1}); continue; end

    tok = regexp(ln, 'Thermal correction to Gibbs Free Energy=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), G_corr = str2double(tok{1}); continue; end

    tok = regexp(ln, 'Sum of electronic and zero-point Energies=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), E0 = str2double(tok{1}); continue; end

    tok = regexp(ln, 'Sum of electronic and thermal Energies=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), U_tot = str2double(tok{1}); continue; end

    tok = regexp(ln, 'Sum of electronic and thermal Enthalpies=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), H_tot = str2double(tok{1}); continue; end

    tok = regexp(ln, 'Sum of electronic and thermal Free Energies=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), G_tot = str2double(tok{1}); continue; end
end

% -------------------------------------------------------------------------
% Assemble output
% -------------------------------------------------------------------------
ha2kJ  = 2625.4996;   % Hartree -> kJ/mol
ha2kcal = 627.5094;   % Hartree -> kcal/mol

en.ZPE_corr  = ZPE_corr;
en.U_corr    = U_corr;
en.H_corr    = H_corr;
en.G_corr    = G_corr;
en.E0        = E0;
en.U         = U_tot;
en.H         = H_tot;
en.G         = G_tot;
en.ZPE_kJ    = ZPE_kJ;
en.T         = T;
en.P         = P_atm;
en.has_thermo = ~isnan(ZPE_corr);

% Entropia: S = (H - G) / T  [Hartree/K] -> [J/(mol·K)]
if en.has_thermo && ~isnan(T) && T > 0
    en.S_JmolK = (H_tot - G_tot) / T * ha2kJ * 1000;
else
    en.S_JmolK = NaN;
end

% Convenience copies in kJ/mol
en.SCF_kJ   = en.SCF    * ha2kJ;
en.G_kJ     = en.G      * ha2kJ;
en.H_kJ     = en.H      * ha2kJ;

en.filename = filename;

% -------------------------------------------------------------------------
% Print summary
% -------------------------------------------------------------------------
fprintf('\n── G16_energy: %s ──\n', filename);
fprintf('  Method : %s\n', en.method);
fprintf('  SCF    : %+.8f  Ha\n', en.SCF);
if en.has_thermo
    fprintf('  ZPE    : %+.8f  Ha   (%.2f kJ/mol)\n', en.ZPE_corr, en.ZPE_kJ);
    fprintf('  E0+ZPE : %+.8f  Ha\n', en.E0);
    fprintf('  H      : %+.8f  Ha\n', en.H);
    fprintf('  G      : %+.8f  Ha\n', en.G);
    fprintf('  S      :  %.4f  J/(mol·K)\n', en.S_JmolK);
    fprintf('  T = %.2f K,  P = %.5f atm\n', en.T, en.P);
end
fprintf('\n');

end
