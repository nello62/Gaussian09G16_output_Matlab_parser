function T = G16_read_all(filename)
% G16_READ_ALL  Collect all relevant data from a Gaussian 16 .out file
%
%   T = G16_READ_ALL(FILENAME) parses a Gaussian 16 output file using the
%   full set of G16_XXX.m toolbox functions and assembles the results into
%   a single struct T, providing a one-call summary of a calculation.
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

    % Mulliken and APT charge values
    T.charge = G16_charges(filename, 'plot', false);

    % SCF etc. data
    T.energy = G16_energy(filename);

    % Optimized structure
    T.structure = G16_structure(filename);

    % Dipole and polarizability data
    T.dipolar = G16_dipole_polar(filename);

    % IR and Raman vibrational modes, if present (freq=Raman)
    T.nmodes = G16_nmodes(filename);

    % Frequency, IR and Raman intensities, IR and Raman spectra (10 cm^-1 FWHM)
    T.spectra = G16_spectra(filename);

    % Info on Gaussian route
    T.route = G16_route(filename);

    % Charge and multiplicity details
    [c, m] = G16_charge_mult(filename);
    T.chargemol.charge = c;
    T.chargemol.mol    = m;

end