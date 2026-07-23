function T = G09_batch_read_all(folder, varargin)
% G09_BATCH_READ_ALL  Runs G09_READ_ALL over every Gaussian 09 output file
%                     in a folder and aggregates the key results into one
%                     summary table.
%
%   T = G09_BATCH_READ_ALL(folder)
%   T = G09_BATCH_READ_ALL(folder, 'Recursive', true)
%   T = G09_BATCH_READ_ALL(folder, 'WriteReports', true)
%   T = G09_BATCH_READ_ALL(folder, 'SaveAs', 'summary.xlsx')
%
%   Scans FOLDER for .log/.out files (case-insensitive), runs
%   G09_READ_ALL (plus G09_ORBITAL_ENERGIES, for the HOMO/LUMO gap) on
%   each, and returns one row per file. A file that fails to parse (e.g.
%   an incomplete/crashed job, or a Gaussian 16 file run through the G09
%   toolbox by mistake) does not stop the batch: its row is filled with
%   NaN and the error message is recorded in .Status instead.
%
%   Optional parameters (Name-Value):
%       'Recursive'    - false (default) | true — also scan subfolders
%       'WriteReports' - false (default) | true — write a
%                        G09_WRITE_REPORT .txt next to each source file
%       'SaveAs'       - '' (default) | path to save the summary table
%                        (.csv or .xlsx, inferred from the extension)
%
%   OUTPUT  table T with one row per file and columns:
%       File, Natoms, SCF_Hartree, E0_Hartree, G_kJmol, mu_tot, mu_units,
%       HOMO_eV, LUMO_eV, Gap_eV, Status
%
%   Example:
%       T = G09_batch_read_all('results/', 'WriteReports', true);
%       T(T.Gap_eV == min(T.Gap_eV), :)   % smallest HOMO-LUMO gap
%
%   See also G09_READ_ALL, G09_WRITE_REPORT, G09_ORBITAL_ENERGIES.

p = inputParser;
addRequired(p,  'folder',       @ischar);
addParameter(p, 'Recursive',    false, @islogical);
addParameter(p, 'WriteReports', false, @islogical);
addParameter(p, 'SaveAs',       '',    @ischar);
parse(p, folder, varargin{:});

recursive     = p.Results.Recursive;
write_reports = p.Results.WriteReports;
save_as       = p.Results.SaveAs;

if ~isfolder(folder)
    error('G09_batch_read_all: folder not found: %s', folder);
end

% -------------------------------------------------------------------------
% Collect .log/.out files (case-insensitive extension match)
% -------------------------------------------------------------------------
if recursive
    d = dir(fullfile(folder, '**', '*'));
else
    d = dir(fullfile(folder, '*'));
end
d = d(~[d.isdir]);

is_target = false(numel(d), 1);
for k = 1:numel(d)
    [~, ~, ext] = fileparts(d(k).name);
    is_target(k) = any(strcmpi(ext, {'.log', '.out'}));
end
d = d(is_target);

if isempty(d)
    error('G09_batch_read_all: no .log/.out files found in %s', folder);
end

Nfiles = numel(d);
File        = cell(Nfiles, 1);
Natoms      = nan(Nfiles, 1);
SCF_Hartree = nan(Nfiles, 1);
E0_Hartree  = nan(Nfiles, 1);
G_kJmol     = nan(Nfiles, 1);
mu_tot      = nan(Nfiles, 1);
mu_units    = cell(Nfiles, 1);
HOMO_eV     = nan(Nfiles, 1);
LUMO_eV     = nan(Nfiles, 1);
Gap_eV      = nan(Nfiles, 1);
Status      = cell(Nfiles, 1);

n_ok = 0;
for k = 1:Nfiles
    fullpath = fullfile(d(k).folder, d(k).name);
    File{k}  = fullpath;
    mu_units{k} = '';

    try
        Tk = G09_read_all(fullpath);
        oe = G09_orbital_energies(fullpath);

        Natoms(k)      = Tk.structure.Natoms;
        SCF_Hartree(k) = Tk.energy.SCF;
        E0_Hartree(k)  = Tk.energy.E0;
        G_kJmol(k)     = Tk.energy.G_kJ;
        mu_tot(k)      = Tk.dipolar.mu_tot;
        mu_units{k}    = Tk.dipolar.mu_units;
        HOMO_eV(k)     = oe.HOMO_eV;
        LUMO_eV(k)     = oe.LUMO_eV;
        Gap_eV(k)      = oe.gap_eV;
        Status{k}      = 'ok';
        n_ok = n_ok + 1;

        if write_reports
            [fdir, fname] = fileparts(fullpath);
            G09_write_report(Tk, fullfile(fdir, [fname, '_report.txt']));
        end
    catch err
        Status{k} = err.message;
    end
end

T = table(File, Natoms, SCF_Hartree, E0_Hartree, G_kJmol, mu_tot, mu_units, ...
    HOMO_eV, LUMO_eV, Gap_eV, Status);

if ~isempty(save_as)
    writetable(T, save_as);
end

fprintf('\n── G09_batch_read_all ──\n');
fprintf('  Folder     : %s\n', folder);
fprintf('  Files found: %d\n', Nfiles);
fprintf('  Succeeded  : %d\n', n_ok);
fprintf('  Failed     : %d\n', Nfiles - n_ok);
if ~isempty(save_as)
    fprintf('  Saved to   : %s\n', save_as);
end
fprintf('\n');

end % G09_batch_read_all
