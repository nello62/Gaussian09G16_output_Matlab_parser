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
%       .nmodes     - vibrational normal modes (see G09_nmodes)
%       .spectra    - IR/Raman spectra (see G09_spectra)
%
%   Example:
%       T = G09_read_all('indaco.log');
%       disp(T.energy)
%       G09_draw_molecule(T.structure);

lines = G09_read_lines(filename);

T.charge    = G09_charges(filename, 'plot', false, 'Lines', lines);
T.energy    = G09_energy(filename, 'Lines', lines);
T.structure = G09_structure(filename, 'Lines', lines);
T.dipolar   = G09_dipole_polar(filename, 'Lines', lines);
T.nmodes    = G09_nmodes(filename, 'Lines', lines);
T.spectra   = G09_spectra(filename, 'Lines', lines);

end % G09_read_all
