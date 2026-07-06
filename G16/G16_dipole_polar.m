function dp = G16_dipole_polar(filename, varargin)
% G16_DIPOLE_POLAR  Estrae momento di dipolo e polarizzabilita' di dipolo
%                   from a Gaussian 16 .out/.log file.
%
%   dp = G16_DIPOLE_POLAR(filename)
%   dp = G16_DIPOLE_POLAR(filename, 'units', 'au')      % default
%   dp = G16_DIPOLE_POLAR(filename, 'units', 'Debye')
%   dp = G16_DIPOLE_POLAR(filename, 'units', 'SI')
%
%   OUTPUT  struct dp with fields:
%   ── Dipole moment ─────────────────────────────────────────────────────
%       .mu_x   .mu_y   .mu_z   Cartesian components
%       .mu_tot                 total magnitude
%       .mu_units               unit string ('au', 'Debye', '10^-30 C·m')
%
%   ── Static polarisability Alpha(0;0) ──────────────────────────────────
%       .alpha_iso              isotropic mean  (1/3)*Tr(alpha)
%       .alpha_aniso            anisotropy
%       .alpha_tensor           [3x3] full tensor (xx,yx,yy,zx,zy,zz)
%       .alpha_units            unit string
%
%   ── Dynamic polarisability Alpha(-w;w) ────────────────────────────────
%       .alpha_dyn              struct array, one entry per laser frequency:
%                                 .lambda_nm   laser wavelength (nm)
%                                 .freq_au     frequency in atomic units
%                                 .iso         isotropic mean
%                                 .aniso       anisotropy
%                                 .tensor      [3x3] polarisability tensor
%       .N_dyn                  number of dynamic (laser) frequencies found
%
%   ── Derivatives with respect to normal modes ──────────────────────────
%       .dmu_dQ    [Nmodes x 3]   dipole derivatives w.r.t. normal modes (dmu_x/dQ, ...)
%                                  (related to IR intensities)
%       .dalpha_dQ [Nmodes x 6]   polarisability derivatives (xx,yx,yy,zx,zy,zz)
%                                  (related to Raman activities)
%       .has_deriv              logical        true if modal derivatives were found
%
%   Note: dmu_dQ and dalpha_dQ are obtained from the same normal mode block;
%         G16 prints them as "RamAct Fr=1, Fr=2" for alpha
%         and as "IR Inten" for |dmu/dQ|². For exact values
%         use G16_spectra which returns correctly scaled IR and Raman activities.
%
%   Example:
%       dp = G16_dipole_polar('V_E00t.out');
%       dp.mu_tot          % modulo dipolo in Debye
%       dp.alpha_iso       % polarizzabilita' isotropica in au
%       dp.alpha_tensor    % full [3x3] polarisability tensor in au

% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename', @ischar);
addParameter(p, 'units',    'au', @ischar);
parse(p, filename, varargin{:});
units = lower(p.Results.units);

if ~isfile(filename)
    error('G16_dipole_polar: file not found: %s', filename);
end
fid  = fopen(filename, 'r');
raw  = fread(fid, '*char')';
fclose(fid);
lines = strsplit(raw, newline);
N = numel(lines);

% -------------------------------------------------------------------------
% Conversions
% -------------------------------------------------------------------------
au2Debye = 2.541746;          % 1 au = 2.541746 D
au2SI_mu = 8.47836e-30;       % 1 au dipolo = 8.478e-30 C·m
au2esu   = 1.48185e-25;       % 1 au alpha = 1.482e-25 cm^3 (= 0.1482 A^3)
au2Ang3  = 0.148185;          % 1 au alpha = 0.1482 Angstrom^3

% -------------------------------------------------------------------------
% Parse: scan file linearly; each quantity is overwritten, so the last occurrence wins
% -------------------------------------------------------------------------
mu_xyz_au   = [NaN NaN NaN];
mu_tot_au   = NaN;
alpha0_iso  = NaN;
alpha0_aniso= NaN;
alpha0_tens = NaN(3,3);   % [xx yx yy zx zy zz] -> 3x3

alpha_dyn_list = struct('lambda_nm', {}, 'freq_au', {}, ...
                        'iso', {}, 'aniso', {}, 'tensor', {});

k = 1;
while k <= N
    ln = lines{k};

    % ── Compact dipole line: "X=  6.1581  Y= -0.0682  Z= -0.0082  Tot=  6.1584"
    if ~isempty(regexp(ln, 'X=.*Y=.*Z=.*Tot=', 'once')) && ...
       ~isempty(strfind(lines{max(1,k-1)}, 'Dipole moment'))
        tok = regexp(ln, 'X=\s*([-\d.]+)\s+Y=\s*([-\d.]+)\s+Z=\s*([-\d.]+)\s+Tot=\s*([\d.]+)', ...
                     'tokens', 'once');
        if ~isempty(tok)
            mu_xyz_au(1) = str2double(tok{1}) / au2Debye;  % D -> au
            mu_xyz_au(2) = str2double(tok{2}) / au2Debye;
            mu_xyz_au(3) = str2double(tok{3}) / au2Debye;
            mu_tot_au    = str2double(tok{4}) / au2Debye;
        end
        k = k+1; continue
    end

    % ── Extended dipole (Fortran format): "x  0.249225D+01  ..."
    % Active only when the previous line is the "Electric dipole moment
    % (input orientation)" header. Gaussian also prints an "(dipole
    % orientation)" variant right after it — a frame rotated so mu aligns
    % with one axis, used internally by the CPHF/polarizability routine —
    % which is intentionally NOT matched here since it does not correspond
    % to the molecular frame used by mol.xyz (G16_structure) elsewhere in
    % this toolbox.
    hdr1 = lines{max(1,k-1)};
    hdr2 = lines{max(1,k-2)};
    is_input_orient = (contains(hdr1, 'Electric dipole moment') && contains(hdr1, 'input orientation')) || ...
                       (contains(hdr2, 'Electric dipole moment') && contains(hdr2, 'input orientation'));
    if is_input_orient
        % reads the Tot, x, y, z values from the following lines
        for off = 0:4
            if k+off > N, break; end
            ln2 = lines{k+off};
            tok_tot = regexp(ln2, '^\s*Tot\s+([-\d.Dd+E]+)', 'tokens', 'once');
            tok_x   = regexp(ln2, '^\s*x\s+([-\d.Dd+E]+)', 'tokens', 'once');
            tok_y   = regexp(ln2, '^\s*y\s+([-\d.Dd+E]+)', 'tokens', 'once');
            tok_z   = regexp(ln2, '^\s*z\s+([-\d.Dd+E]+)', 'tokens', 'once');
            if ~isempty(tok_tot), mu_tot_au  = fortran2double(tok_tot{1}); end
            if ~isempty(tok_x),   mu_xyz_au(1) = fortran2double(tok_x{1}); end
            if ~isempty(tok_y),   mu_xyz_au(2) = fortran2double(tok_y{1}); end
            if ~isempty(tok_z),   mu_xyz_au(3) = fortran2double(tok_z{1}); end
        end
    end

    % ── Alpha(0;0): static block
    if ~isempty(regexp(ln, '^\s*Alpha\(0;0\):', 'once'))
        % Extract iso, aniso, and tensor elements from the following lines
        alpha_tmp = parse_alpha_block(lines, k+1, N);
        if ~isempty(alpha_tmp)
            alpha0_iso   = alpha_tmp.iso;
            alpha0_aniso = alpha_tmp.aniso;
            alpha0_tens  = alpha_tmp.tensor;
        end
        k = k + 10; continue
    end

    % ── Alpha(-w;w): dynamic block
    tok_dyn = regexp(ln, 'Alpha\(-w;w\)\s+w=\s*([\d.]+)nm', 'tokens', 'once');
    if isempty(tok_dyn)
        tok_dyn = regexp(ln, 'Alpha\(-w;w\)\s+w=\s*([\d.]+)\s*nm', 'tokens', 'once');
    end
    if ~isempty(tok_dyn)
        lam_nm = str2double(tok_dyn{1});
        alpha_tmp = parse_alpha_block(lines, k+1, N);
        if ~isempty(alpha_tmp)
            nd = numel(alpha_dyn_list) + 1;
            alpha_dyn_list(nd).lambda_nm = lam_nm;
            alpha_dyn_list(nd).freq_au   = 45.5640 / lam_nm;  % nm -> au (appross.)
            alpha_dyn_list(nd).iso       = alpha_tmp.iso;
            alpha_dyn_list(nd).aniso     = alpha_tmp.aniso;
            alpha_dyn_list(nd).tensor    = alpha_tmp.tensor;
        end
        k = k + 10; continue
    end

    k = k + 1;
end

% -------------------------------------------------------------------------
% Apply unit conversion factors
% -------------------------------------------------------------------------
switch units
    case 'au'
        ufac_mu    = 1;
        ufac_alpha = 1;
        ulbl_mu    = 'au';
        ulbl_alpha = 'au (bohr^3)';
    case 'debye'
        ufac_mu    = au2Debye;
        ufac_alpha = 1;          % alpha resta in au
        ulbl_mu    = 'Debye';
        ulbl_alpha = 'au (bohr^3)';
    case 'si'
        ufac_mu    = au2SI_mu;
        ufac_alpha = au2esu * 1e6;   % au -> 10^-24 esu (cm^3) come in G16
        ulbl_mu    = '10^-30 C·m';
        ulbl_alpha = '10^-24 esu (cm^3)';
    otherwise
        error('G16_dipole_polar: units must be ''au'', ''Debye'' or ''SI''.');
end

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
dp.mu_x     = mu_xyz_au(1) * ufac_mu;
dp.mu_y     = mu_xyz_au(2) * ufac_mu;
dp.mu_z     = mu_xyz_au(3) * ufac_mu;
dp.mu_tot   = mu_tot_au    * ufac_mu;
dp.mu_units = ulbl_mu;

dp.alpha_iso    = alpha0_iso   * ufac_alpha;
dp.alpha_aniso  = alpha0_aniso * ufac_alpha;
dp.alpha_tensor = alpha0_tens  * ufac_alpha;
dp.alpha_units  = ulbl_alpha;

% Apply unit conversion to dynamic entries
for nd = 1:numel(alpha_dyn_list)
    alpha_dyn_list(nd).iso    = alpha_dyn_list(nd).iso    * ufac_alpha;
    alpha_dyn_list(nd).aniso  = alpha_dyn_list(nd).aniso  * ufac_alpha;
    alpha_dyn_list(nd).tensor = alpha_dyn_list(nd).tensor * ufac_alpha;
end
dp.alpha_dyn  = alpha_dyn_list;
dp.N_dyn      = numel(alpha_dyn_list);

dp.has_deriv  = false;   % derivate modali via G16_spectra (IR/Raman)
dp.filename   = filename;

% -------------------------------------------------------------------------
% Print summary
% -------------------------------------------------------------------------
fprintf('\n── G16_dipole_polar: %s ──\n', filename);
fprintf('  Dipolo  μ = (%.4f, %.4f, %.4f)  |μ| = %.4f %s\n', ...
    dp.mu_x, dp.mu_y, dp.mu_z, dp.mu_tot, dp.mu_units);
fprintf('  α iso   = %.3f %s\n', dp.alpha_iso,   dp.alpha_units);
fprintf('  α aniso = %.3f %s\n', dp.alpha_aniso, dp.alpha_units);
fprintf('  α tensor [au]:\n');
fprintf('    %10.3f  %10.3f  %10.3f\n', alpha0_tens(1,1), alpha0_tens(1,2), alpha0_tens(1,3));
fprintf('    %10.3f  %10.3f  %10.3f\n', alpha0_tens(2,1), alpha0_tens(2,2), alpha0_tens(2,3));
fprintf('    %10.3f  %10.3f  %10.3f\n', alpha0_tens(3,1), alpha0_tens(3,2), alpha0_tens(3,3));
if dp.N_dyn > 0
    fprintf('  Alpha(-w;w):\n');
    for nd = 1:dp.N_dyn
        fprintf('    λ=%.1f nm  iso=%.3f  aniso=%.3f\n', ...
            dp.alpha_dyn(nd).lambda_nm, dp.alpha_dyn(nd).iso, dp.alpha_dyn(nd).aniso);
    end
end
fprintf('\n');

end  % G16_dipole_polar


% =========================================================================
%  Local function: read iso/aniso/xx/yx/yy/zx/zy/zz block
%  Returns a struct with .iso, .aniso, .tensor [3x3], or [] if parsing fails
% =========================================================================
function out = parse_alpha_block(lines, k_start, N)
out = [];
iso = NaN; aniso = NaN;
xx = NaN; yx = NaN; yy = NaN; zx = NaN; zy = NaN; zz = NaN;

for k = k_start : min(k_start+10, N)
    ln = lines{k};
    % Line format: "   label   value_au   value_esu   value_SI"
    % The au value is the first number after the label
    tok = regexp(ln, '^\s*(\w+)\s+([-\d.Dd+E]+)', 'tokens', 'once');
    if isempty(tok), continue; end
    label = lower(tok{1});
    val   = fortran2double(tok{2});
    switch label
        case 'iso',   iso   = val;
        case 'aniso', aniso = val;
        case 'xx',    xx    = val;
        case 'yx',    yx    = val;
        case 'yy',    yy    = val;
        case 'zx',    zx    = val;
        case 'zy',    zy    = val;
        case 'zz',    zz    = val;
    end
end

if isnan(xx), return; end

% Assemble the symmetric [3x3] tensor (alpha is symmetric by definition)
out.iso    = iso;
out.aniso  = aniso;
out.tensor = [xx, yx, zx; ...
              yx, yy, zy; ...
              zx, zy, zz];
end


% =========================================================================
%  Local function: convert Fortran D notation -> MATLAB double
% =========================================================================
function v = fortran2double(s)
v = str2double(strrep(strrep(s, 'D', 'e'), 'd', 'e'));
end
