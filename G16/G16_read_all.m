function T = G16_read_all(filename)
% G16_READ_ALL  Collect all relevant data from a Gaussian 16 .out file
%
%   T = G16_READ_ALL(FILENAME) parses a Gaussian 16 output file using the
%   full set of G16_XXX.m toolbox functions and assembles the results into
%   a single struct T, providing a one-call summary of a calculation. The
%   file is read from disk once and the parsed lines are reused by every
%   sub-function (via their 'Lines' parameter), instead of each one
%   re-reading and re-splitting the file independently.
%
%   INPUT:
%     filename - path to the Gaussian 16 .out (or .log) file
%
%   OUTPUT:
%     T - struct with fields:
%       .charge          - Mulliken and APT atomic charges
%                           (see G16_charges; plotting disabled)
%       .energy          - SCF energy and related energetics
%                           (see G16_energy)
%       .structure       - optimized molecular structure
%                           (see G16_structure)
%       .dipolar         - dipole moment and polarizability data
%                           (see G16_dipole_polar)
%       .nmodes          - IR/Raman vibrational normal modes (see
%                           G16_nmodes); field entirely omitted, with a
%                           non-blocking warning, if the file has no
%                           frequency calculation (e.g. a TD-DFT-only
%                           single point job with no 'freq' keyword) --
%                           check with isfield(T, 'nmodes')
%       .spectra         - frequencies, IR/Raman intensities, and
%                           simulated IR/Raman spectra (10 cm^-1 FWHM)
%                           (see G16_spectra); omitted under the same
%                           condition as .nmodes above
%       .route           - Gaussian route section details
%                           (see G16_route)
%       .chargemol.charge - total molecular charge
%       .chargemol.mol    - spin multiplicity
%                           (see G16_charge_mult)
%
%   EXAMPLE:
%     T = G16_read_all('zeatin.out');
%     disp(T.energy)
%     G09_draw_molecule(T.structure);
%
%   S.Trusso IPCF-CNR Messina 2026
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

    if ~isfile(filename)
        error('G16_read_all: file not found: %s', filename);
    end
    fid  = fopen(filename, 'r');
    raw  = fread(fid, '*char')';
    fclose(fid);
    lines = strsplit(raw, newline);
    G16_check_gaussian_match(lines, filename);

    % Mulliken and APT charge values
    T.charge = G16_charges(filename, 'plot', false, 'Lines', lines);

    % SCF etc. data
    T.energy = G16_energy(filename, 'Lines', lines);

    % Optimized structure
    T.structure = G16_structure(filename, 'Lines', lines);

    % Dipole and polarizability data
    T.dipolar = G16_dipole_polar(filename, 'Lines', lines);

    % IR and Raman vibrational modes, if present (requires 'freq' in the
    % route section, e.g. absent for a TD-DFT-only single point job) --
    % caught here rather than left to crash the whole call, since
    % G16_write_report already expects T.nmodes to be optional (isfield check)
    try
        T.nmodes = G16_nmodes(filename, 'Lines', lines);
    catch ME
        warning('G16_read_all:noNmodes', ...
            'No vibrational normal modes found in %s (%s); T.nmodes omitted.', ...
            filename, ME.message);
    end

    % Frequency, IR and Raman intensities, IR and Raman spectra (10 cm^-1 FWHM)
    % -- same optional-section handling as T.nmodes above
    try
        T.spectra = G16_spectra(filename, 'Lines', lines);
    catch ME
        warning('G16_read_all:noSpectra', ...
            'No IR/Raman spectra found in %s (%s); T.spectra omitted.', ...
            filename, ME.message);
    end

    % Info on Gaussian route
    T.route = G16_route(filename, 'Lines', lines);

    % Charge and multiplicity details
    [c, m] = G16_charge_mult(filename, 'Lines', lines);
    T.chargemol.charge = c;
    T.chargemol.mol    = m;

end
