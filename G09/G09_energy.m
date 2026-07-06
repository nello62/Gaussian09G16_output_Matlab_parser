function en = G09_energy(filename, varargin)
% G09_ENERGY  Extracts SCF energy and thermochemistry from a Gaussian 09 output file.
%
%   en = G09_ENERGY(filename)
%   en = G09_ENERGY(filename, 'step', 'last')
%   en = G09_ENERGY(filename, 'step', N)
%
%   OUTPUT  struct en with fields (all in Hartree unless stated):
%       .SCF          SCF energy (Hartree)
%       .method       method string, e.g. 'RB3LYP'
%       .ZPE_corr     zero-point correction (Hartree)
%       .U_corr       thermal correction to energy (Hartree)
%       .H_corr       thermal correction to enthalpy (Hartree)
%       .G_corr       thermal correction to Gibbs free energy (Hartree)
%       .E0           SCF + ZPE (Hartree)
%       .H            SCF + H_corr (Hartree)
%       .G            SCF + G_corr (Hartree)
%       .S_JmolK      entropy S (J mol-1 K-1)
%       .ZPE_kJ       zero-point energy (kJ/mol)
%       .T            temperature (K)
%       .P            pressure (atm)
%       .has_thermo   logical
%       .G_kJ         G in kJ/mol
%       .H_kJ         H in kJ/mol
%       .filename     char

% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename', @ischar);
addParameter(p, 'step',     'last', @(x) ischar(x) || isnumeric(x));
parse(p, filename, varargin{:});
step_req = p.Results.step;

lines = G09_read_lines(filename);

% -------------------------------------------------------------------------
% Collect all SCF Done lines
% -------------------------------------------------------------------------
scf_vals    = [];
scf_methods = {};

for k = 1:numel(lines)
    ln = lines{k};
    tok = regexp(ln, 'SCF Done:\s+E\((\S+)\)\s*=\s*([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok)
        scf_methods{end+1} = tok{1}; %#ok<AGROW>
        scf_vals(end+1)    = str2double(tok{2}); %#ok<AGROW>
    end
end

if isempty(scf_vals)
    error('G09_energy: no "SCF Done" line found in %s', filename);
end

n_scf = numel(scf_vals);
if ischar(step_req)
    si = n_scf - 1 + (strcmp(lower(step_req),'last') - strcmp(lower(step_req),'last'));
    % simpler:
    if strcmp(lower(step_req), 'last'),  si = n_scf;
    elseif strcmp(lower(step_req), 'first'), si = 1;
    else, error('G09_energy: step must be ''first'', ''last'', or an integer.');
    end
else
    si = round(step_req);
    if si < 1 || si > n_scf
        error('G09_energy: step %d out of range [1, %d].', si, n_scf);
    end
end

SCF    = scf_vals(si);
method = scf_methods{si};

% -------------------------------------------------------------------------
% Thermochemistry (last occurrence — identical format to G16)
% -------------------------------------------------------------------------
ZPE_corr = NaN; U_corr = NaN; H_corr = NaN; G_corr = NaN;
E0 = NaN; H_tot = NaN; G_tot = NaN;
ZPE_kJ = NaN; T = NaN; P_atm = NaN;

for k = 1:numel(lines)
    ln = lines{k};

    tok = regexp(ln, 'Temperature\s+([\d.]+)\s+Kelvin', 'tokens', 'once');
    if ~isempty(tok), T = str2double(tok{1}); end

    tok = regexp(ln, 'Pressure\s+([\d.]+)', 'tokens', 'once');
    if ~isempty(tok), P_atm = str2double(tok{1}); end

    tok = regexp(ln, 'Zero-point vibrational energy\s+([\d.]+)\s+\(Joules', 'tokens', 'once');
    if ~isempty(tok), ZPE_kJ = str2double(tok{1}) / 1000; end

    tok = regexp(ln, 'Zero-point correction=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), ZPE_corr = str2double(tok{1}); end

    tok = regexp(ln, 'Thermal correction to Energy=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), U_corr = str2double(tok{1}); end

    tok = regexp(ln, 'Thermal correction to Enthalpy=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), H_corr = str2double(tok{1}); end

    tok = regexp(ln, 'Thermal correction to Gibbs Free Energy=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), G_corr = str2double(tok{1}); end

    tok = regexp(ln, 'Sum of electronic and zero-point Energies=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), E0 = str2double(tok{1}); end

    tok = regexp(ln, 'Sum of electronic and thermal Enthalpies=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), H_tot = str2double(tok{1}); end

    tok = regexp(ln, 'Sum of electronic and thermal Free Energies=\s+([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), G_tot = str2double(tok{1}); end
end

ha2kJ = 2625.4996;
has_thermo = ~isnan(ZPE_corr);
S_JmolK = NaN;
if has_thermo && ~isnan(T) && T > 0 && ~isnan(H_tot) && ~isnan(G_tot)
    S_JmolK = (H_tot - G_tot) / T * ha2kJ * 1000;
end

en.SCF        = SCF;
en.method     = method;
en.ZPE_corr   = ZPE_corr;
en.U_corr     = U_corr;
en.H_corr     = H_corr;
en.G_corr     = G_corr;
en.E0         = E0;
en.H          = H_tot;
en.G          = G_tot;
en.S_JmolK    = S_JmolK;
en.ZPE_kJ     = ZPE_kJ;
en.T          = T;
en.P          = P_atm;
en.has_thermo = has_thermo;
en.G_kJ       = G_tot * ha2kJ;
en.H_kJ       = H_tot * ha2kJ;
en.filename   = filename;

fprintf('\n── G09_energy: %s ──\n', filename);
fprintf('  Method : %s\n', method);
fprintf('  SCF    : %+.8f  Ha\n', SCF);
if has_thermo
    fprintf('  ZPE    : %+.8f  Ha  (%.2f kJ/mol)\n', ZPE_corr, ZPE_kJ);
    fprintf('  H      : %+.8f  Ha\n', H_tot);
    fprintf('  G      : %+.8f  Ha\n', G_tot);
    fprintf('  S      :  %.4f  J/(mol·K)\n', S_JmolK);
    fprintf('  T = %.2f K,  P = %.5f atm\n', T, P_atm);
end
fprintf('\n');

end  % G09_energy
