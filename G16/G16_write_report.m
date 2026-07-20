function outfile = G16_write_report(T, outfile)
% G16_WRITE_REPORT  Writes a human-readable text report from a G16_read_all struct.
%
%   G16_WRITE_REPORT(T) writes a formatted summary of every field in T
%   (the struct returned by G16_READ_ALL) to a .txt file, named after the
%   source Gaussian file (e.g. 'zeatin.out' -> 'zeatin_report.txt') in the
%   current folder.
%
%   G16_WRITE_REPORT(T, OUTFILE) writes to OUTFILE instead.
%
%   OUTFILE = G16_WRITE_REPORT(...) also returns the path written to.
%
%   Sections included (only if the corresponding field is present in T):
%   route, charge/multiplicity, molecular structure, energetics, dipole
%   moment and polarisability, atomic charges, vibrational modes, and a
%   summary of the simulated IR/Raman spectra (the full continuum arrays
%   are not dumped — use T.spectra.x/.IR_cont/.Raman_cont directly for that).
%
%   Example:
%       T = G16_read_all('zeatin.out');
%       G16_write_report(T);
%       G16_write_report(T, 'zeatin_summary.txt');
%
%   See also G16_READ_ALL.

if nargin < 2 || isempty(outfile)
    src = '';
    if isfield(T, 'structure') && isfield(T.structure, 'filename')
        src = T.structure.filename;
    end
    if isempty(src)
        outfile = 'G16_report.txt';
    else
        [~, name] = fileparts(src);
        outfile = [name '_report.txt'];
    end
end

fid = fopen(outfile, 'w');
if fid == -1
    error('G16_write_report:cannotOpen', 'Could not open %s for writing.', outfile);
end
cleaner = onCleanup(@() fclose(fid));

fprintf(fid, '================================================================\n');
fprintf(fid, '  Gaussian 16 Calculation Report\n');
fprintf(fid, '================================================================\n');
if isfield(T, 'structure') && isfield(T.structure, 'filename')
    fprintf(fid, 'Source file : %s\n', T.structure.filename);
end
fprintf(fid, 'Generated   : %s\n\n', datestr(now));

% -------------------------------------------------------------------
if isfield(T, 'route')
    fprintf(fid, '--- Route section ---------------------------------------------\n');
    fprintf(fid, '%s\n\n', T.route);
end

% -------------------------------------------------------------------
if isfield(T, 'chargemol')
    fprintf(fid, '--- Charge / multiplicity ---------------------------------------\n');
    fprintf(fid, 'Total charge       : %d\n', T.chargemol.charge);
    fprintf(fid, 'Spin multiplicity  : %d\n\n', T.chargemol.mol);
end

% -------------------------------------------------------------------
if isfield(T, 'structure')
    s = T.structure;
    fprintf(fid, '--- Molecular structure -----------------------------------------\n');
    fprintf(fid, 'Atoms       : %d\n', s.Natoms);
    if isfield(s, 'orientation'), fprintf(fid, 'Orientation : %s\n', s.orientation); end
    if isfield(s, 'step'),        fprintf(fid, 'Step        : %d\n', s.step); end
    fprintf(fid, '\n%-6s %-4s %12s %12s %12s\n', 'Idx', 'Sym', 'X (A)', 'Y (A)', 'Z (A)');
    for i = 1:s.Natoms
        fprintf(fid, '%-6d %-4s %12.6f %12.6f %12.6f\n', ...
            i, s.symbols{i}, s.xyz(i,1), s.xyz(i,2), s.xyz(i,3));
    end
    fprintf(fid, '\n');
end

% -------------------------------------------------------------------
if isfield(T, 'energy')
    e = T.energy;
    fprintf(fid, '--- Energetics ----------------------------------------------------\n');
    fprintf(fid, 'Method          : %s\n', e.method);
    fprintf(fid, 'SCF energy      : %.8f Hartree\n', e.SCF);
    if e.has_thermo
        fprintf(fid, 'ZPE correction  : %.8f Hartree  (%.3f kJ/mol)\n', e.ZPE_corr, e.ZPE_kJ);
        fprintf(fid, 'Thermal U corr. : %.8f Hartree\n', e.U_corr);
        fprintf(fid, 'Thermal H corr. : %.8f Hartree\n', e.H_corr);
        fprintf(fid, 'Thermal G corr. : %.8f Hartree\n', e.G_corr);
        fprintf(fid, 'E0 (SCF+ZPE)    : %.8f Hartree\n', e.E0);
        if isfield(e, 'U')
            u_val = e.U;
        else
            u_val = e.SCF + e.U_corr;   % defensive: compute if missing
        end
        fprintf(fid, 'U               : %.8f Hartree\n', u_val);
        fprintf(fid, 'H               : %.8f Hartree\n', e.H);
        fprintf(fid, 'G               : %.8f Hartree\n', e.G);
        fprintf(fid, 'T, P            : %.2f K, %.4f atm\n', e.T, e.P);
    else
        fprintf(fid, '(no thermochemistry data — opt/single-point job without freq)\n');
    end
    fprintf(fid, '\n');
end

% -------------------------------------------------------------------
if isfield(T, 'dipolar')
    d = T.dipolar;
    fprintf(fid, '--- Dipole moment and polarisability -------------------------------\n');
    fprintf(fid, 'Dipole (mu_x, mu_y, mu_z) : %.6f  %.6f  %.6f  [%s]\n', ...
        d.mu_x, d.mu_y, d.mu_z, d.mu_units);
    fprintf(fid, 'Dipole magnitude          : %.6f %s\n', d.mu_tot, d.mu_units);
    fprintf(fid, 'Alpha isotropic           : %.6f %s\n', d.alpha_iso, d.alpha_units);
    fprintf(fid, 'Alpha anisotropy          : %.6f %s\n', d.alpha_aniso, d.alpha_units);
    if isfield(d, 'N_dyn') && d.N_dyn > 0
        fprintf(fid, '\nDynamic polarisability Alpha(-w;w):\n');
        fprintf(fid, '%-14s %-14s %14s %14s\n', 'Lambda (nm)', 'Freq (au)', 'Iso', 'Aniso');
        for k = 1:d.N_dyn
            ad = d.alpha_dyn(k);
            fprintf(fid, '%-14.2f %-14.6f %14.6f %14.6f\n', ...
                ad.lambda_nm, ad.freq_au, ad.iso, ad.aniso);
        end
    end
    fprintf(fid, '\n');
end

% -------------------------------------------------------------------
if isfield(T, 'charge')
    c = T.charge;
    fprintf(fid, '--- Atomic charges (%s) --------------------------------------------\n', c.type);
    fprintf(fid, '%-6s %-4s %14s\n', 'Idx', 'Sym', 'Charge (e)');
    for i = 1:c.Natoms
        fprintf(fid, '%-6d %-4s %14.6f\n', i, c.symbols{i}, c.charges(i));
    end
    fprintf(fid, 'Sum of charges : %.6f\n', c.sum_q);
    if isfield(c, 'dipole') && ~isempty(c.dipole)
        fprintf(fid, 'Dipole (from charges overlay) : %.6f  %.6f  %.6f  Debye\n', c.dipole);
    end
    fprintf(fid, '\n');
end

% -------------------------------------------------------------------
if isfield(T, 'nmodes')
    nm = T.nmodes;
    fprintf(fid, '--- Vibrational normal modes ---------------------------------------\n');
    fprintf(fid, 'Number of modes : %d\n\n', nm.Nmodes);
    if nm.has_Raman
        fprintf(fid, '%-6s %12s %10s %12s %10s %8s\n', ...
            'Mode', 'Freq(cm-1)', 'IR', 'Raman', 'RedMass', 'Sym');
    else
        fprintf(fid, '%-6s %12s %10s %10s %8s\n', ...
            'Mode', 'Freq(cm-1)', 'IR', 'RedMass', 'Sym');
    end
    for k = 1:nm.Nmodes
        sym = '';
        if numel(nm.symmetry) >= k, sym = nm.symmetry{k}; end
        if nm.has_Raman
            fprintf(fid, '%-6d %12.2f %10.2f %12.2f %10.4f %8s\n', ...
                k, nm.freq(k), nm.IR(k), nm.Raman(k), nm.redmass(k), sym);
        else
            fprintf(fid, '%-6d %12.2f %10.2f %10.4f %8s\n', ...
                k, nm.freq(k), nm.IR(k), nm.redmass(k), sym);
        end
    end
    fprintf(fid, '\n');
end

% -------------------------------------------------------------------
if isfield(T, 'spectra')
    sp = T.spectra;
    fprintf(fid, '--- Simulated IR/Raman spectra --------------------------------------\n');
    fprintf(fid, 'FWHM used     : %.2f cm^-1\n', sp.FWHM);
    fprintf(fid, 'Grid range    : %.1f - %.1f cm^-1 (%d points)\n', ...
        sp.x(1), sp.x(end), numel(sp.x));
    fprintf(fid, 'Raman present : %s\n', mat2str(sp.has_Raman));
    fprintf(fid, '(full continuum arrays not dumped here — see T.spectra.x / .IR_cont / .Raman_cont)\n\n');
end

fprintf('G16_write_report: report written to %s\n', outfile);

end % G16_write_report
