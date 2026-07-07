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
%       .nmodes          - IR/Raman vibrational normal modes, if present
%                           (see G16_nmodes)
%       .spectra         - frequencies, IR/Raman intensities, and
%                           simulated IR/Raman spectra (10 cm^-1 FWHM)
%                           (see G16_spectra)
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

    if ~isfile(filename)
        error('G16_read_all: file not found: %s', filename);
    end
    fid  = fopen(filename, 'r');
    raw  = fread(fid, '*char')';
    fclose(fid);
    lines = strsplit(raw, newline);

    % Mulliken and APT charge values
    T.charge = G16_charges(filename, 'plot', false, 'Lines', lines);

    % SCF etc. data
    T.energy = G16_energy(filename, 'Lines', lines);

    % Optimized structure
    T.structure = G16_structure(filename, 'Lines', lines);

    % Dipole and polarizability data
    T.dipolar = G16_dipole_polar(filename, 'Lines', lines);

    % IR and Raman vibrational modes, if present (freq=Raman)
    T.nmodes = G16_nmodes(filename, 'Lines', lines);

    % Frequency, IR and Raman intensities, IR and Raman spectra (10 cm^-1 FWHM)
    T.spectra = G16_spectra(filename, 'Lines', lines);

    % Info on Gaussian route
    T.route = G16_route(filename, 'Lines', lines);

    % Charge and multiplicity details
    [c, m] = G16_charge_mult(filename, 'Lines', lines);
    T.chargemol.charge = c;
    T.chargemol.mol    = m;

end
