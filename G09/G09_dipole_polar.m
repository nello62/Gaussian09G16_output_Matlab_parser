function dp = G09_dipole_polar(filename, varargin)
% G09_DIPOLE_POLAR  Extracts dipole moment and polarisability from a Gaussian 09 file.
%
%   dp = G09_DIPOLE_POLAR(filename)
%   dp = G09_DIPOLE_POLAR(filename, 'units', 'Debye')
%
%   DIFFERENCES FROM G16:
%     - G09 prints "Exact polarizability: xx yx yy zx zy zz" on a SINGLE line
%       (6 values, upper triangle in row order).
%       G16 prints formatted Alpha(0;0) / Alpha(-w;w) blocks.
%     - G09 may also print "Approx polarizability: ..." (CPHF approximation);
%       this function reads only the "Exact" value.
%     - G09 has NO dynamic alpha entries (no laser frequency dependence)
%       unless explicitly requested with Polar=OptRot or similar keywords.
%     - G09 also prints "Diagonal vibrational polarizability:" (3 values).
%
%   Optional parameters:
%       'units'  - 'au' (default) | 'Debye' (for mu only) | 'SI'
%       'Lines'  - pre-read cell array of file lines (from G09_READ_LINES),
%                  to skip re-reading the file when it has already been
%                  read elsewhere (e.g. G09_READ_ALL). Default {} (read
%                  the file normally).
%
%   OUTPUT  struct dp with fields:
%       .mu_x .mu_y .mu_z   dipole components
%       .mu_tot             dipole magnitude
%       .mu_units           unit string
%       .alpha_iso          isotropic polarisability = (1/3)Tr(alpha)
%       .alpha_aniso        anisotropy
%       .alpha_tensor       [3x3] tensor (au)
%       .alpha_units        unit string
%       .alpha_vib          [1x3] diagonal vibrational polarisability (au)
%       .has_alpha_vib      logical
%       .alpha_approx       [3x3] CPHF approximate tensor (au), NaN if absent
%       .filename           char

% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename', @ischar);
addParameter(p, 'units',    'au', @ischar);
addParameter(p, 'Lines',    {},   @iscell);
parse(p, filename, varargin{:});
units = lower(p.Results.units);

lines = p.Results.Lines;
if isempty(lines)
    if ~isfile(filename)
        error('G09_dipole_polar: file not found: %s', filename);
    end
    lines = G09_read_lines(filename);
end
N     = numel(lines);

% Conversion factors
au2Debye = 2.541746;
switch units
    case 'au',    ufac_mu = 1;          ulbl_mu = 'au';
    case 'debye', ufac_mu = au2Debye;   ulbl_mu = 'Debye';
    case 'si',    ufac_mu = 8.47836e-30; ulbl_mu = '10^{-30} C·m';
    otherwise,    error('G09_dipole_polar: units must be ''au'', ''Debye'', or ''SI''.');
end

% -------------------------------------------------------------------------
% Parse
% -------------------------------------------------------------------------
mu_xyz_D      = [NaN NaN NaN];
mu_tot_D      = NaN;
alpha_exact   = NaN(3,3);
alpha_approx  = NaN(3,3);
alpha_vib     = NaN(1,3);

for k = 1:N
    ln = lines{k};

    % ── Dipole moment (Debye)
    % Format: "    X=    -5.8857    Y=     0.1068    Z=     0.0065  Tot=    5.8867"
    if k > 1 && ~isempty(strfind(lines{k-1}, 'Dipole moment'))
        tok = regexp(ln, 'X=\s*([-\d.]+)\s+Y=\s*([-\d.]+)\s+Z=\s*([-\d.]+)\s+Tot=\s*([\d.]+)', ...
                     'tokens', 'once');
        if ~isempty(tok)
            mu_xyz_D = [str2double(tok{1}), str2double(tok{2}), str2double(tok{3})];
            mu_tot_D = str2double(tok{4});
        end
    end

    % ── Exact polarizability (6 values: xx yx yy zx zy zz)
    % Format: "  Exact polarizability: 465.591  -7.755 211.415  -0.007   0.004 100.077"
    if ~isempty(strfind(ln, 'Exact polarizability:'))
        vals = sscanf(strrep(ln, 'Exact polarizability:', ''), '%f');
        if numel(vals) >= 6
            % upper triangle: (1,1)=xx (2,1)=yx (2,2)=yy (3,1)=zx (3,2)=zy (3,3)=zz
            alpha_exact = build_tensor_from_upper(vals(1:6));
        end
    end

    % ── Approx polarizability (optional, same format)
    if ~isempty(strfind(ln, 'Approx polarizability:'))
        vals = sscanf(strrep(ln, 'Approx polarizability:', ''), '%f');
        if numel(vals) >= 6
            alpha_approx = build_tensor_from_upper(vals(1:6));
        end
    end

    % ── Diagonal vibrational polarizability (3 values on next line)
    if ~isempty(strfind(ln, 'Diagonal vibrational polarizability:'))
        if k+1 <= N
            vals = sscanf(lines{k+1}, '%f');
            if numel(vals) >= 3
                alpha_vib = vals(1:3)';
            end
        end
    end
end

% -------------------------------------------------------------------------
% Derived quantities
% -------------------------------------------------------------------------
if ~isnan(alpha_exact(1,1))
    trace_a     = alpha_exact(1,1) + alpha_exact(2,2) + alpha_exact(3,3);
    alpha_iso   = trace_a / 3;
    alpha_aniso = sqrt(0.5 * ( ...
        (alpha_exact(1,1)-alpha_exact(2,2))^2 + ...
        (alpha_exact(2,2)-alpha_exact(3,3))^2 + ...
        (alpha_exact(3,3)-alpha_exact(1,1))^2 + ...
        6*(alpha_exact(1,2)^2 + alpha_exact(1,3)^2 + alpha_exact(2,3)^2)));
else
    alpha_iso   = NaN;
    alpha_aniso = NaN;
end

mu_xyz = mu_xyz_D / au2Debye * ufac_mu;
mu_tot = mu_tot_D / au2Debye * ufac_mu;

dp.mu_x         = mu_xyz(1);
dp.mu_y         = mu_xyz(2);
dp.mu_z         = mu_xyz(3);
dp.mu_tot       = mu_tot;
dp.mu_units     = ulbl_mu;
dp.alpha_iso    = alpha_iso;
dp.alpha_aniso  = alpha_aniso;
dp.alpha_tensor = alpha_exact;
dp.alpha_units  = 'au (bohr^3)';
dp.alpha_approx = alpha_approx;
dp.alpha_vib    = alpha_vib;
dp.has_alpha_vib = ~any(isnan(alpha_vib));
dp.filename     = filename;

fprintf('\n── G09_dipole_polar: %s ──\n', filename);
fprintf('  μ = (%.4f, %.4f, %.4f)  |μ| = %.4f %s\n', ...
    dp.mu_x, dp.mu_y, dp.mu_z, dp.mu_tot, ulbl_mu);
fprintf('  α_iso   = %.3f au\n', alpha_iso);
fprintf('  α_aniso = %.3f au\n', alpha_aniso);
fprintf('  α tensor (au):\n');
fprintf('    %10.3f  %10.3f  %10.3f\n', alpha_exact(1,1), alpha_exact(1,2), alpha_exact(1,3));
fprintf('    %10.3f  %10.3f  %10.3f\n', alpha_exact(2,1), alpha_exact(2,2), alpha_exact(2,3));
fprintf('    %10.3f  %10.3f  %10.3f\n', alpha_exact(3,1), alpha_exact(3,2), alpha_exact(3,3));
if dp.has_alpha_vib
    fprintf('  α_vib (diag) = [%.4f  %.4f  %.4f] au\n', alpha_vib(1), alpha_vib(2), alpha_vib(3));
end
fprintf('\n');

end  % G09_dipole_polar


% =========================================================================
%  Local: build symmetric 3x3 from upper triangle (row order)
%  Input order: xx yx yy zx zy zz
% =========================================================================
function T = build_tensor_from_upper(v)
% v = [xx, yx, yy, zx, zy, zz]
T = [v(1), v(2), v(4); ...
     v(2), v(3), v(5); ...
     v(4), v(5), v(6)];
end
