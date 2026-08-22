function result = AU_convert(value, quantity, direction, varargin)
% AU_CONVERT  Converts physical quantities between atomic units (au) and SI
%             (or other common units used in computational chemistry).
%
%   result = AU_CONVERT(value, quantity, direction)
%   result = AU_CONVERT(value, quantity, direction, 'target', unit)
%
%   All conversion factors are based on CODATA 2022 (NIST).
%
% ─────────────────────────────────────────────────────────────────────────
%  INPUTS
% ─────────────────────────────────────────────────────────────────────────
%   value     : numeric scalar or array to convert
%   quantity  : string — physical quantity (see table below)
%   direction : 'au2si'  — from atomic units to SI (or target unit)
%               'si2au'  — from SI (or target unit) to atomic units
%
%   Optional Name-Value:
%   'target'  : target unit (only needed when multiple SI-side units exist)
%               Default target for each quantity is listed below.
%
% ─────────────────────────────────────────────────────────────────────────
%  SUPPORTED QUANTITIES  (quantity string, case-insensitive)
% ─────────────────────────────────────────────────────────────────────────
%
%  ENERGY
%   'energy'        au (Hartree) <-> J           (default)
%   'energy'        au (Hartree) <-> kJ/mol      (target: 'kJ/mol')
%   'energy'        au (Hartree) <-> kcal/mol    (target: 'kcal/mol')
%   'energy'        au (Hartree) <-> eV          (target: 'eV')
%   'energy'        au (Hartree) <-> cm-1        (target: 'cm-1')
%
%  LENGTH
%   'length'        au (Bohr)   <-> m            (default)
%   'length'        au (Bohr)   <-> Angstrom     (target: 'Angstrom')
%   'length'        au (Bohr)   <-> pm           (target: 'pm')
%
%  DIPOLE MOMENT
%   'dipole'        au (e·a0)   <-> C·m          (default)
%   'dipole'        au (e·a0)   <-> Debye        (target: 'Debye')
%
%  POLARISABILITY
%   'polar'         au          <-> C²·m²·J⁻¹   (default, SI)
%   'polar'         au          <-> Å³           (target: 'Angstrom3')
%   'polar'         au          <-> cm³ (esu)    (target: 'esu')
%
%  FIRST HYPERPOLARISABILITY
%   'beta'          au          <-> C³·m³·J⁻²   (default, SI ×10⁻⁵³)
%   'beta'          au          <-> esu           (target: 'esu', ×10⁻³³)
%
%  FORCE
%   'force'         au (Eh/a0)  <-> N            (default)
%   'force'         au (Eh/a0)  <-> nN           (target: 'nN')
%
%  FORCE CONSTANT
%   'frcconst'      au (Eh/a0²) <-> N/m          (default)
%   'frcconst'      au (Eh/a0²) <-> mDyne/Å      (target: 'mDyne/A')
%
%  FREQUENCY
%   'frequency'     au (Eh/ħ)   <-> rad/s        (default)
%   'frequency'     au (Eh/ħ)   <-> Hz           (target: 'Hz')
%   'frequency'     au (Eh/ħ)   <-> cm-1         (target: 'cm-1')
%   'frequency'     au (Eh/ħ)   <-> THz          (target: 'THz')
%
%  TIME
%   'time'          au (ħ/Eh)   <-> s            (default)
%   'time'          au (ħ/Eh)   <-> fs           (target: 'fs')
%   'time'          au (ħ/Eh)   <-> as           (target: 'as')
%
%  PRESSURE
%   'pressure'      au (Eh/a0³) <-> Pa           (default)
%   'pressure'      au (Eh/a0³) <-> GPa          (target: 'GPa')
%   'pressure'      au (Eh/a0³) <-> atm          (target: 'atm')
%
%  ELECTRIC FIELD
%   'efield'        au (Eh/ea0) <-> V/m          (default)
%   'efield'        au (Eh/ea0) <-> GV/m         (target: 'GV/m')
%   'efield'        au (Eh/ea0) <-> V/Angstrom   (target: 'V/A')
%
%  MASS
%   'mass'          au (me)     <-> kg           (default)
%   'mass'          au (me)     <-> Da (AMU)     (target: 'Da')
%
%  CHARGE
%   'charge'        au (e)      <-> C            (default)
%
%  MAGNETIC MOMENT
%   'magmom'        au (μB)     <-> J/T          (default)
%
%  TEMPERATURE  (via thermal energy kBT = Eh)
%   'temp'          K           <-> Eh           (direction 'si2au': K->Eh)
%
% ─────────────────────────────────────────────────────────────────────────
%  EXAMPLES
% ─────────────────────────────────────────────────────────────────────────
%   AU_convert(-875.932, 'energy',   'au2si', 'target', 'kJ/mol')
%   AU_convert(6.3347,   'dipole',   'si2au', 'target', 'Debye')
%   AU_convert(259.03,   'polar',    'au2si', 'target', 'Angstrom3')
%   AU_convert(1582.8,   'frequency','si2au', 'target', 'cm-1')
%   AU_convert(0.76,     'length',   'si2au', 'target', 'Angstrom')
%   AU_convert(mol.xyz,  'length',   'au2si', 'target', 'Angstrom')  % array
%
%   % List all supported quantities:
%   AU_convert([], 'help', '')
%
% ─────────────────────────────────────────────────────────────────────────

% ── Help mode ─────────────────────────────────────────────────────────────
if strcmpi(quantity, 'help')
    print_help();
    result = [];
    return
end

% ── Parse optional target ─────────────────────────────────────────────────
p = inputParser;
addRequired(p, 'value');
addRequired(p, 'quantity', @ischar);
addRequired(p, 'direction', @ischar);
addParameter(p, 'target', '', @ischar);
parse(p, value, quantity, direction, varargin{:});
target = lower(strtrim(p.Results.target));
dir    = lower(strtrim(direction));

if ~ismember(dir, {'au2si','si2au'})
    error('AU_convert: direction must be ''au2si'' or ''si2au''.');
end

% ── CODATA 2022 constants ─────────────────────────────────────────────────
e    = 1.602176634e-19;       % C            (exact)
me   = 9.1093837139e-31;      % kg
hbar = 1.054571817e-34;       % J·s          (exact)
h    = 6.62607015e-34;        % J·s          (exact)
c    = 299792458.0;           % m/s          (exact)
kB   = 1.380649e-23;          % J/K          (exact)
NA   = 6.02214076e23;         % mol-1        (exact)
a0   = 5.29177210544e-11;     % m            (Bohr radius, CODATA 2022)
Eh   = 4.3597447222060e-18;   % J            (Hartree, CODATA 2022)
eps0 = 8.8541878188e-12;      % F/m

% ── Dispatch ──────────────────────────────────────────────────────────────
q = lower(strtrim(quantity));

switch q

  % ── ENERGY ──────────────────────────────────────────────────────────────
  case 'energy'
    switch target
        case {'kj/mol','kjmol','kj'}
            factor = Eh * NA / 1000;  % 2625.499639 kJ/mol
            unit_si = 'kJ/mol';
        case {'kcal/mol','kcalmol','kcal'}
            factor = Eh * NA / 4184;  % 627.509474 kcal/mol
            unit_si = 'kcal/mol';
        case {'ev','electronvolt'}
            factor = Eh / e;          % 27.211386 eV
            unit_si = 'eV';
        case {'cm-1','cm^-1','wavenumber','wn'}
            factor = Eh / (h * c * 100);  % 219474.631 cm-1
            unit_si = 'cm-1';
        otherwise  % default: J
            factor = Eh;
            unit_si = 'J';
    end

  % ── LENGTH ──────────────────────────────────────────────────────────────
  case 'length'
    switch target
        case {'angstrom','a','ang','å'}
            factor = a0 * 1e10;       % 0.529177 Å
            unit_si = 'Å';
        case {'pm','picometer'}
            factor = a0 * 1e12;       % 52.9177 pm
            unit_si = 'pm';
        otherwise  % default: m
            factor = a0;
            unit_si = 'm';
    end

  % ── DIPOLE MOMENT ───────────────────────────────────────────────────────
  case {'dipole','dipolemoment','mu'}
    switch target
        case {'debye','d'}
            factor = (e * a0) / 3.33564095e-30;  % 2.541747 D
            unit_si = 'Debye';
        otherwise  % default: C·m
            factor = e * a0;          % 8.478354e-30 C·m
            unit_si = 'C·m';
    end

  % ── POLARISABILITY ──────────────────────────────────────────────────────
  case {'polar','polarisability','polarizability','alpha'}
    % 1 au = e^2 a0^2 / Eh  (SI: C^2 m^2 J^-1)
    % also = a0^3 / (4 pi eps0) in m^3 volume units
    alpha_si = e^2 * a0^2 / Eh;   % 1.648777e-41 C^2 m^2 J^-1
    switch target
        case {'angstrom3','a3','ang3','å3'}
            % volume polarisability: alpha/(4 pi eps0) in Å^3
            factor = alpha_si / (4*pi*eps0) * 1e30;  % 0.148185 Å^3
            unit_si = 'Å³';
        case {'esu','cm3'}
            % CGS volume: alpha/(4 pi eps0) in cm^3  x 1e6
            factor = alpha_si / (4*pi*eps0) * 1e6;   % 1.48185e-25 cm^3
            unit_si = 'cm³ (esu)';
        otherwise  % default: C^2 m^2 J^-1 (SI)
            factor = alpha_si;
            unit_si = 'C²·m²·J⁻¹';
    end

  % ── FIRST HYPERPOLARISABILITY ────────────────────────────────────────────
  case {'beta','hyperpolar','firsthyperpolar'}
    % 1 au = e^3 a0^3 / Eh^2
    beta_si = e^3 * a0^3 / Eh^2;  % 3.2064e-53 C^3 m^3 J^-2
    switch target
        case {'esu'}
            % 1 au = 8.6392e-33 esu (statC^3 cm^3 erg^-2)
            factor = 8.6392e-33;
            unit_si = '×10⁻³³ esu';
        otherwise  % default: C^3 m^3 J^-2 (SI) expressed as ×10^-53
            factor = beta_si;
            unit_si = 'C³·m³·J⁻²';
    end

  % ── SECOND HYPERPOLARISABILITY ───────────────────────────────────────────
  case {'gamma','secondhyperpolar'}
    % 1 au = e^4 a0^4 / Eh^3
    gamma_si = e^4 * a0^4 / Eh^3;
    factor  = gamma_si;
    unit_si = 'C⁴·m⁴·J⁻³';

  % ── FORCE ────────────────────────────────────────────────────────────────
  case 'force'
    F_au = Eh / a0;               % 8.238724e-8 N
    switch target
        case {'nn','nanonewton'}
            factor = F_au * 1e9;
            unit_si = 'nN';
        otherwise
            factor = F_au;
            unit_si = 'N';
    end

  % ── FORCE CONSTANT ───────────────────────────────────────────────────────
  case {'frcconst','forceconstant','frc'}
    % 1 au = Eh/a0^2
    k_au = Eh / a0^2;             % 1556.893 N/m
    switch target
        case {'mdyne/a','mdyne/ang','mdyne/angstrom','mdyne'}
            % 1 N/m = 1e-2 mDyne/Å → 1 au = k_au * 1e-2 mDyne/Å
            factor = k_au * 1e-2;  % 15.5689 mDyne/Å
            unit_si = 'mDyne/Å';
        otherwise
            factor = k_au;
            unit_si = 'N/m';
    end

  % ── FREQUENCY ────────────────────────────────────────────────────────────
  case {'frequency','freq','omega'}
    omega_au = Eh / hbar;         % 4.134137e16 rad/s
    switch target
        case {'hz','hertz'}
            factor = omega_au / (2*pi);
            unit_si = 'Hz';
        case {'thz','terahertz'}
            factor = omega_au / (2*pi) / 1e12;
            unit_si = 'THz';
        case {'cm-1','cm^-1','wavenumber','wn'}
            factor = omega_au / (2*pi*c*100);  % 219474.6 cm-1
            unit_si = 'cm⁻¹';
        otherwise  % rad/s
            factor = omega_au;
            unit_si = 'rad/s';
    end

  % ── TIME ─────────────────────────────────────────────────────────────────
  case 'time'
    t_au = hbar / Eh;             % 2.418884e-17 s
    switch target
        case {'fs','femtosecond'}
            factor = t_au * 1e15;
            unit_si = 'fs';
        case {'as','attosecond'}
            factor = t_au * 1e18;
            unit_si = 'as';
        otherwise
            factor = t_au;
            unit_si = 's';
    end

  % ── PRESSURE ─────────────────────────────────────────────────────────────
  case 'pressure'
    P_au = Eh / a0^3;             % 2.942102e13 Pa
    switch target
        case {'gpa','gigapascal'}
            factor = P_au / 1e9;
            unit_si = 'GPa';
        case {'atm','atmosphere'}
            factor = P_au / 101325;
            unit_si = 'atm';
        case {'bar'}
            factor = P_au / 1e5;
            unit_si = 'bar';
        otherwise
            factor = P_au;
            unit_si = 'Pa';
    end

  % ── ELECTRIC FIELD ───────────────────────────────────────────────────────
  case {'efield','electricfield'}
    E_au = Eh / (e * a0);         % 5.14221e11 V/m
    switch target
        case {'gv/m','gvm'}
            factor = E_au / 1e9;
            unit_si = 'GV/m';
        case {'v/a','v/angstrom','v/ang'}
            factor = E_au * 1e-10;  % V/Å
            unit_si = 'V/Å';
        case {'mv/cm','mvcm'}
            factor = E_au / 1e6;
            unit_si = 'MV/cm';
        otherwise
            factor = E_au;
            unit_si = 'V/m';
    end

  % ── MASS ─────────────────────────────────────────────────────────────────
  case 'mass'
    switch target
        case {'da','amu','u','dalton'}
            % 1 au (me) = me/u = 5.48580e-4 Da
            factor = me / 1.66053906892e-27;
            unit_si = 'Da (AMU)';
        otherwise
            factor = me;
            unit_si = 'kg';
    end

  % ── CHARGE ───────────────────────────────────────────────────────────────
  case 'charge'
    factor  = e;
    unit_si = 'C';

  % ── MAGNETIC MOMENT ──────────────────────────────────────────────────────
  case {'magmom','magneticmoment'}
    % 1 au = Bohr magneton = e*hbar/(2*me)
    mu_B   = e * hbar / (2 * me);   % 9.2740100657e-24 J/T
    factor  = mu_B;
    unit_si = 'J/T';

  % ── TEMPERATURE (kBT) ────────────────────────────────────────────────────
  case {'temp','temperature'}
    % Direction: si2au means K -> Eh (thermal energy kBT/Eh)
    % Direction: au2si means Eh -> K
    factor  = kB / Eh;             % 3.16681e-6 Eh/K
    unit_si = 'K (via kBT)';

  otherwise
    error(['AU_convert: unknown quantity ''%s''.\n' ...
           'Call AU_convert([], ''help'', '''') for the full list.'], quantity);
end

% ── Apply conversion ──────────────────────────────────────────────────────
if strcmp(dir, 'au2si')
    result = value .* factor;
else
    result = value ./ factor;
end

% ── Print result if single value ──────────────────────────────────────────
if isscalar(value)
    if strcmp(dir, 'au2si')
        fprintf('%g au  →  %g %s\n', value, result, unit_si);
    else
        fprintf('%g %s  →  %g au\n', value, unit_si, result);
    end
end

end  % AU_convert


% =========================================================================
%  Print help table
% =========================================================================
function print_help()
fprintf('\n');
fprintf('══════════════════════════════════════════════════════════════════\n');
fprintf('  AU_CONVERT — Atomic Units <-> SI Conversion Toolbox\n');
fprintf('  CODATA 2022 constants  |  Sebastiano Trusso, CNR-IPCF Messina\n');
fprintf('══════════════════════════════════════════════════════════════════\n\n');
fprintf('  Usage: result = AU_convert(value, quantity, direction, ''target'', unit)\n\n');

rows = {
  'QUANTITY',    'au UNIT',    'SI/target',           'FACTOR (au->SI)';
  '─────────',  '──────────', '──────────────────',  '────────────────';
  'energy',      'Hartree',    'J',                   '4.359745e-18';
  '',            '',           'kJ/mol',              '2625.4996';
  '',            '',           'kcal/mol',            '627.5095';
  '',            '',           'eV',                  '27.211386';
  '',            '',           'cm-1',                '219474.631';
  'length',      'Bohr (a0)', 'm',                   '5.291772e-11';
  '',            '',           'Angstrom',            '0.529177';
  '',            '',           'pm',                  '52.9177';
  'dipole',      'e·a0',      'C·m',                 '8.478354e-30';
  '',            '',           'Debye',               '2.541747';
  'polar',       'au',         'C²·m²·J⁻¹',          '1.648777e-41';
  '',            '',           'Angstrom3',           '0.148185';
  '',            '',           'esu (cm³)',           '1.48185e-25';
  'beta',        'au',         'C³·m³·J⁻²',          '3.206361e-53';
  '',            '',           'esu',                 '8.6392e-33';
  'gamma',       'au',         'C⁴·m⁴·J⁻³',          '6.235380e-65';
  'force',       'Eh/a0',     'N',                   '8.238724e-8';
  '',            '',           'nN',                  '82.38724';
  'frcconst',    'Eh/a0²',    'N/m',                 '1556.893';
  '',            '',           'mDyne/Å',             '15.5689';
  'frequency',   'Eh/ħ',      'rad/s',               '4.134137e16';
  '',            '',           'Hz',                  '6.579684e15';
  '',            '',           'THz',                 '6579.684';
  '',            '',           'cm-1',                '219474.631';
  'time',        'ħ/Eh',      's',                   '2.418884e-17';
  '',            '',           'fs',                  '0.024189';
  '',            '',           'as',                  '24.1888';
  'pressure',    'Eh/a0³',    'Pa',                  '2.942102e13';
  '',            '',           'GPa',                 '29421.0';
  '',            '',           'atm',                 '2.904430e8';
  'efield',      'Eh/(e·a0)', 'V/m',                 '5.142207e11';
  '',            '',           'GV/m',                '514.221';
  '',            '',           'V/Å',                 '51.4221';
  'mass',        'me',         'kg',                  '9.109384e-31';
  '',            '',           'Da (AMU)',             '5.485799e-4';
  'charge',      'e',          'C',                   '1.602177e-19';
  'magmom',      'μB',         'J/T',                 '9.274010e-24';
  'temp',        'Eh',         'K  (via kBT)',        '3.157750e5';
};

for i = 1:size(rows, 1)
    fprintf('  %-14s %-12s %-22s %s\n', rows{i,:});
end

fprintf('\nExamples:\n');
fprintf('  AU_convert(-875.932, ''energy'',    ''au2si'', ''target'', ''kJ/mol'')\n');
fprintf('  AU_convert(6.3347,   ''dipole'',    ''si2au'', ''target'', ''Debye'')\n');
fprintf('  AU_convert(259.03,   ''polar'',     ''au2si'', ''target'', ''Angstrom3'')\n');
fprintf('  AU_convert(1582.8,   ''frequency'', ''si2au'', ''target'', ''cm-1'')\n');
fprintf('  AU_convert(mol.xyz,  ''length'',    ''au2si'', ''target'', ''Angstrom'')\n');
fprintf('  AU_convert([],       ''help'',      '''')\n\n');
end