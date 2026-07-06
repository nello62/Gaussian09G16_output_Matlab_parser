T=G09_read_all(filename)

T.charge=G09_charges(filename);
T.energy=G09_energy(filename);
T.structure=G09_structure(filename);
T.dipolar=G09_dipole_polar(filename);
T.nmodes=G09_nmodes(filename);
T.spectra=G09_spectra(filename);
T