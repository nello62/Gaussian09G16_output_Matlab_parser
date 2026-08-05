function T = G09_available_data(filename)
% G09_AVAILABLE_DATA  Lists which toolbox quantities/functions are
%                     expected to be extractable from a Gaussian 09
%                     output file, based on the route section keywords
%                     actually present.
%
%   T = G09_AVAILABLE_DATA(filename)
%
%   Reads the route section (see G09_ROUTE) and checks it against the
%   keyword requirements documented in the toolbox manual ("Which .in
%   keywords generate which .out section"), so you can tell up front
%   which data will actually be present in the file -- e.g. an opt-only
%   job (no 'freq') has no vibrational normal modes, no IR/Raman
%   spectra, and no thermochemistry (ZPE/H/G/S), so calling
%   G09_nmodes/G09_spectra on it raises an error and G09_energy's
%   *_corr/E0/U/H/G/T/P/S fields come back NaN -- none of that is a bug,
%   it is simply not in the file. This function lets you check that in
%   one call before running the extraction functions, rather than
%   discovering it one error/NaN at a time.
%
%   OUTPUT  table T with columns:
%       Function    string   toolbox function/quantity name
%       Available   logical  true if the route contains the keyword(s)
%                            needed to generate that data
%       Requires    string   human-readable keyword requirement
%       Notes       string   short clarifying note
%
%   Example:
%       G09_available_data('indaco.log');
%       T = G09_available_data('molecule.out');
%       T(~T.Available, :)   % see what's NOT available in this file
%
%   Note: this predicts availability from the route keywords alone; it
%   does not itself read the data sections, so an incomplete/crashed job
%   can still show a keyword as present without the corresponding output
%   actually having been printed (the extraction function will still
%   raise its own clear error in that case).
%
%   See also G09_ROUTE, G09_READ_ALL.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

route = G09_route(filename);

tokens = regexp(lower(route), '\S+', 'match');

has_freq        = local_has_kw(tokens, 'freq');
has_opt         = local_has_kw(tokens, 'opt');
has_polar       = local_has_kw(tokens, 'polar');
has_freq_raman  = any(cellfun(@(t) (startsWith(t, 'freq=') || startsWith(t, 'freq(')) && contains(t, 'raman'), tokens));
has_cphf_rdfreq = any(cellfun(@(t) startsWith(t, 'cphf=') && contains(t, 'rdfreq'), tokens));
has_nbo         = any(cellfun(@(t) startsWith(t, 'pop=') && contains(t, 'nbo'), tokens));

Function = {};
Available = [];
Requires = {};
Notes = {};

    function add(name, avail, req, note)
        Function{end+1}  = name;
        Available(end+1) = avail;
        Requires{end+1}  = req;
        Notes{end+1}     = note;
    end

add('structure',                       true,           '(default)',                      'nosymm/nosym suppresses Standard orientation only');
add('energy (SCF)',                    true,           '(default)',                      '');
add('energy (thermochemistry: ZPE/H/G/S)', has_freq,   'freq',                           '');
add('charges (Mulliken)',              true,           '(default)',                      '');
add('charges (APT)',                   has_freq,       'freq',                           'printed automatically with any freq job');
add('dipole_polar (dipole moment)',    true,           '(default)',                      '');
add('dipole_polar (static Alpha)',     has_polar,      'polar',                          '');
add('dipole_polar (dynamic Alpha)',    has_polar && has_cphf_rdfreq, 'polar + cphf=rdfreq', 'also needs an explicit frequency list in the .in file');
add('nmodes',                          has_freq,       'freq',                           '');
add('spectra (IR)',                    has_freq,       'freq',                           '');
add('spectra (Raman)',                 has_freq_raman, 'freq=raman',                     '');
add('orbital_energies',                true,           '(default)',                      'full eigenvalue list printed unconditionally');
add('convergence',                     has_opt,        'opt',                            '');
add('charge_mult',                     true,           '(default)',                      '');
add('route',                           true,            '(default)',                     '');
add('get_bond_length',                 true,           '(default, derived)',             'purely geometric, not a real Gaussian section');
add('nbo_bonds',                       has_nbo,        'pop=nbo',                        '');
add('gaussian_version',                true,           '(default)',                      '');
add('restart',                         true,           '(any completed step)',           '');

T = table(Function(:), logical(Available(:)), Requires(:), Notes(:), ...
    'VariableNames', {'Function', 'Available', 'Requires', 'Notes'});

fprintf('\n── G09_available_data: %s ──\n', filename);
fprintf('  Route: %s\n\n', route);
name_w = max(cellfun(@length, Function));
req_w  = max(cellfun(@length, Requires));
for i = 1:height(T)
    if T.Available(i)
        mark = 'YES';
    else
        mark = 'no ';
    end
    fprintf('  %s  %-*s  needs: %-*s  %s\n', mark, name_w, Function{i}, req_w, Requires{i}, Notes{i});
end
fprintf('\n');

end % G09_available_data


% =========================================================================
%  Local function: does any token match keyword "name" as itself,
%  "name=...", or "name(...)"?
% =========================================================================
function tf = local_has_kw(tokens, name)
tf = any(cellfun(@(t) strcmp(t, name) || startsWith(t, [name '=']) || startsWith(t, [name '(']), tokens));
end
