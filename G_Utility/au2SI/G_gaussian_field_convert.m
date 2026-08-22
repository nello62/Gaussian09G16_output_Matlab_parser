function out = G_gaussian_field_convert(value, direction, varargin)
% G_GAUSSIAN_FIELD_CONVERT  Converts between Gaussian's Field=X+N route
%                           keyword integer and a physical electric-field
%                           strength (V/Angstrom, V/m, V/cm, or a.u.).
%
%   out = G_GAUSSIAN_FIELD_CONVERT(N, 'g2phys')
%   out = G_GAUSSIAN_FIELD_CONVERT(N, 'g2phys', 'Unit', 'V/Ang')
%   out = G_GAUSSIAN_FIELD_CONVERT(value, 'phys2g', 'Unit', 'au')
%
%   Gaussian's Field route keyword (e.g. Field=X+10, Field=X-5) specifies
%   a static electric field whose magnitude, in atomic units, is the
%   keyword's integer N times 0.0001 -- i.e. N*0.0001 a.u. Verified
%   against Gaussian's own keyword reference (gaussian.com/field/):
%   "N*0.0001 specifies the magnitude of the field in atomic units in the
%   first format" (the M+-N multipole-field format, which is what a plain
%   dipole field like Field=X+10 uses).
%
%   direction:
%       'g2phys' - value is the Gaussian keyword integer N (e.g. the "10"
%                  in Field=X+10); returns the field strength in 'Unit'
%                  (default 'au').
%       'phys2g' - value is a field strength in 'Unit' (default 'au');
%                  returns the Gaussian keyword integer N, rounded to the
%                  nearest integer -- Gaussian's keyword can only express
%                  multiples of 0.0001 a.u., so a value that is not
%                  already such a multiple is necessarily rounded. A
%                  warning is printed if that rounding changes the value
%                  by more than 0.5%.
%
%   Optional parameters:
%       'Unit'  - 'au' | 'V/m' | 'V/cm' | 'V/Ang' (default 'au')
%       'Print' - print a one-line summary (default true)
%
%   value may be a numeric array; the conversion is applied elementwise.
%
%   Example:
%       G_gaussian_field_convert(10, 'g2phys')                   % 0.0010 (au)
%       G_gaussian_field_convert(10, 'g2phys', 'Unit', 'V/Ang')  % 0.5142 V/Ang
%       G_gaussian_field_convert(0.0025, 'phys2g', 'Unit', 'au') % 25
%       G_gaussian_field_convert(0.5, 'phys2g', 'Unit', 'V/Ang') % N for a 0.5 V/Ang field
%
%   See also AU_CONVERT.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

p = inputParser;
addRequired(p,  'value', @isnumeric);
addRequired(p,  'direction', @(x) any(strcmpi(x, {'g2phys', 'phys2g'})));
addParameter(p, 'Unit',  'au',  @(x) any(strcmpi(x, {'au', 'V/m', 'V/cm', 'V/Ang'})));
addParameter(p, 'Print', true,  @islogical);
parse(p, value, direction, varargin{:});

direction = lower(p.Results.direction);
unit      = p.Results.Unit;
do_print  = p.Results.Print;

% -------------------------------------------------------------------------
% CODATA 2022 constants (same values used throughout this toolbox, e.g.
% AU_convert's own 'efield' quantity) -- the atomic unit of electric field
% strength is computed from them the same way, for exact numerical
% consistency rather than a separately hardcoded rounded literal.
% -------------------------------------------------------------------------
e_C   = 1.602176634e-19;      % C  (exact)
a0_m  = 5.29177210544e-11;    % m  (Bohr radius)
Eh_J  = 4.3597447222060e-18;  % J  (Hartree)
au2Vm = Eh_J / (e_C * a0_m);  % V/m per atomic unit of electric field

switch lower(unit)
    case 'au'
        unit_per_au = 1;
        unit_label  = 'au';
    case 'v/m'
        unit_per_au = au2Vm;
        unit_label  = 'V/m';
    case 'v/cm'
        unit_per_au = au2Vm / 100;         % 1 m = 100 cm
        unit_label  = 'V/cm';
    case 'v/ang'
        unit_per_au = au2Vm * 1e-10;       % 1 Angstrom = 1e-10 m
        unit_label  = 'V/Ang';
end

% -------------------------------------------------------------------------
% Convert
% -------------------------------------------------------------------------
switch direction
    case 'g2phys'
        field_au = value * 1e-4;           % Gaussian's N*0.0001 convention
        out = field_au * unit_per_au;
        if do_print
            fprintf('G_gaussian_field_convert: Field keyword N=%s  ->  %s %s\n', ...
                mat2str(value), mat2str(out, 6), unit_label);
        end

    case 'phys2g'
        field_au = value / unit_per_au;
        N_raw = field_au / 1e-4;
        out = round(N_raw);
        rel_err = abs(out - N_raw) ./ max(abs(N_raw), eps);
        if do_print
            fprintf('G_gaussian_field_convert: %s %s  ->  Field keyword N=%s\n', ...
                mat2str(value), unit_label, mat2str(out));
        end
        if any(rel_err(:) > 0.005)
            warning('G_gaussian_field_convert:rounded', ...
                ['The requested field strength is not an exact multiple of 0.0001 a.u. ', ...
                 '(Gaussian''s Field keyword granularity); rounded to the nearest integer N, ', ...
                 'a relative change of up to %.2f%%.'], 100*max(rel_err(:)));
        end
end

end % G_gaussian_field_convert
