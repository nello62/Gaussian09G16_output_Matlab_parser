function T = G09_read_all(filename)
% G09_READ_ALL  Collect all relevant data from a Gaussian 09 output file.
%
%   T = G09_READ_ALL(filename) parses a Gaussian 09 output file using the
%   full set of G09_XXX.m toolbox functions and assembles the results into
%   a single struct T, providing a one-call summary of a calculation. The
%   file is read from disk once and the parsed lines are reused by every
%   sub-function (via their 'Lines' parameter), instead of each one
%   re-reading and re-splitting the file independently.
%
%   OUTPUT  struct T with fields:
%       .charge     - Mulliken and APT atomic charges (see G09_charges; plotting disabled)
%       .energy     - SCF energy and thermochemistry (see G09_energy)
%       .structure  - molecular geometry (see G09_structure)
%       .dipolar    - dipole moment and polarisability (see G09_dipole_polar)
%       .nmodes     - vibrational normal modes (see G09_nmodes); field
%                     entirely omitted, with a non-blocking warning, if
%                     the file has no frequency calculation (e.g. a
%                     TD-DFT-only single point job with no 'freq'
%                     keyword) -- check with isfield(T, 'nmodes')
%       .spectra    - IR/Raman spectra (see G09_spectra); omitted under
%                     the same condition as .nmodes above
%
%   Example:
%       T = G09_read_all('indaco.log');
%       disp(T.energy)
%       G09_draw_molecule(T.structure);
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

lines = G09_read_lines(filename);

T.charge    = G09_charges(filename, 'plot', false, 'Lines', lines);
T.energy    = G09_energy(filename, 'Lines', lines);
T.structure = G09_structure(filename, 'Lines', lines);
T.dipolar   = G09_dipole_polar(filename, 'Lines', lines);

% IR and Raman vibrational modes, if present (requires 'freq' in the
% route section, e.g. absent for a TD-DFT-only single point job) --
% caught here rather than left to crash the whole call, since
% G09_write_report already expects T.nmodes to be optional (isfield check)
try
    T.nmodes = G09_nmodes(filename, 'Lines', lines);
catch ME
    warning('G09_read_all:noNmodes', ...
        'No vibrational normal modes found in %s (%s); T.nmodes omitted.', ...
        filename, ME.message);
end

% Frequency, IR and Raman intensities, IR and Raman spectra -- same
% optional-section handling as T.nmodes above
try
    T.spectra = G09_spectra(filename, 'Lines', lines);
catch ME
    warning('G09_read_all:noSpectra', ...
        'No IR/Raman spectra found in %s (%s); T.spectra omitted.', ...
        filename, ME.message);
end

end % G09_read_all
