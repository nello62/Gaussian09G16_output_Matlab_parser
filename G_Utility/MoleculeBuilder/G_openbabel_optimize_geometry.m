function xyzOut = G_openbabel_optimize_geometry(symbols, xyz, bondList, varargin)
%G_OPENBABEL_OPTIMIZE_GEOMETRY  Geometry optimization via Open Babel's
%   force-field minimizer -- a more accurate clean-up than
%   G_ROUGH_OPTIMIZE_GEOMETRY's simple geometric penalty, for a hand-built
%   molecule (e.g. from G_MOLECULE_BUILDER) before sending it to Gaussian
%   for a real quantum-chemical optimization.
%
%   xyzOut = G_OPENBABEL_OPTIMIZE_GEOMETRY(symbols, xyz, bondList)
%   xyzOut = G_OPENBABEL_OPTIMIZE_GEOMETRY(..., Name, Value, ...)
%
%   Requires Open Babel (https://openbabel.org) installed and reachable --
%   the only function in this toolbox with an external, non-MATLAB
%   dependency (everything else is self-contained by design). If you do
%   not have it:
%       conda create -n obabel_env -c conda-forge openbabel -y
%   (or any other install method -- Homebrew, a system package manager,
%   a manual build). This function auto-detects a few common install
%   locations (see LOCATE_OBABEL below); pass 'ObabelPath' explicitly if
%   none of them match.
%
%   Input:
%       symbols  - {Natoms x 1 cell} atomic symbols
%       xyz      - [Natoms x 3 double] starting Cartesian coordinates (Angstrom)
%       bondList - [Nbonds x 2] or [Nbonds x 3] Atom1/Atom2(/Order) pairs
%                  (same convention as G16_DRAW_MOLECULE's 'BondList'; a
%                  missing 3rd column defaults every bond to order 1)
%
%   Optional parameters (Name-Value):
%       'ForceField' - 'MMFF94' (default) | 'MMFF94s' | 'UFF' | 'GAFF' |
%                      'Ghemical' -- any force field Open Babel's own
%                      `obabel -L forcefields` lists. MMFF94 covers
%                      organic molecules well; UFF covers a much wider
%                      element range (most of the periodic table) at
%                      lower accuracy -- use it for anything MMFF94
%                      rejects (it errors on elements it has no
%                      parameters for, e.g. most metals).
%       'Steps'      - maximum minimization steps (default: 2500, Open
%                      Babel's own default)
%       'ObabelPath' - full path to the `obabel` executable (default: '',
%                      auto-detected)
%
%   How it works: the molecule is written to a temporary MDL Molfile
%   (.mol, V2000 -- carries explicit bonds/bond orders, so Open Babel
%   optimizes the exact connectivity BONDLIST specifies, rather than
%   re-perceiving bonds from geometry the way a bare .xyz input would
%   force it to), Open Babel is invoked as
%   `obabel in.mol -O out.xyz --minimize --ff <ForceField> --noh --steps <Steps>`,
%   and the resulting geometry is read back from out.xyz. '--noh' is
%   always passed: Open Babel otherwise adds explicit hydrogens to
%   complete any open valence before minimizing, which would silently
%   change the atom count on return -- this function always returns
%   exactly as many atoms, in the same order, as it was given (an error
%   is raised instead, in the unexpected case that this invariant is
%   ever violated).
%
%   A molecule with fewer than 2 atoms is returned unchanged (nothing for
%   Open Babel to optimize).
%
%   See also G_ROUGH_OPTIMIZE_GEOMETRY, G_MOLECULE_BUILDER, G16_DRAW_MOLECULE.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

    p = inputParser;
    addRequired(p,  'symbols');
    addRequired(p,  'xyz',      @isnumeric);
    addRequired(p,  'bondList', @isnumeric);
    addParameter(p, 'ForceField', 'MMFF94', @ischar);
    addParameter(p, 'Steps',      2500,     @isnumeric);
    addParameter(p, 'ObabelPath', '',       @ischar);
    parse(p, symbols, xyz, bondList, varargin{:});
    opt = p.Results;

    symbols = cellstr(symbols(:));
    N = numel(symbols);
    if size(xyz, 1) ~= N
        error('G_openbabel_optimize_geometry:sizeMismatch', ...
            'Number of symbols (%d) does not match number of xyz rows (%d).', N, size(xyz, 1));
    end
    if N < 2
        xyzOut = xyz;
        return
    end

    obabelPath = opt.ObabelPath;
    if isempty(obabelPath)
        obabelPath = locateObabel();
    end

    tmpDir = tempname;
    mkdir(tmpDir);
    cleanupObj = onCleanup(@() rmdir(tmpDir, 's'));  % deletes the temp dir on return, including on error
    molFile = fullfile(tmpDir, 'in.mol');
    xyzFile = fullfile(tmpDir, 'out.xyz');

    writeMolfileV2000(molFile, symbols, xyz, bondList);

    cmd = sprintf('%s %s -O %s --minimize --ff %s --noh --steps %d', ...
        quoteForShell(obabelPath), quoteForShell(molFile), quoteForShell(xyzFile), ...
        opt.ForceField, round(opt.Steps));
    [status, cmdout] = system(cmd);
    % Open Babel can exit 0 even when it fails outright (confirmed: e.g.
    % an element the chosen force field has no parameters for -- "Could
    % not setup force field" -- still exits 0 and writes an EMPTY output
    % file), so a zero exit code and an existing file are not enough on
    % their own; its own "N molecule(s) converted" summary line is the
    % reliable signal, checked here instead.
    nConverted = regexp(cmdout, '(\d+)\s+molecules?\s+converted', 'tokens', 'once');
    convertedOk = ~isempty(nConverted) && str2double(nConverted{1}) >= 1;
    if status ~= 0 || ~isfile(xyzFile) || ~convertedOk || dir(xyzFile).bytes == 0
        error('G_openbabel_optimize_geometry:obabelFailed', ...
            'Open Babel failed (exit code %d): %s', status, strtrim(cmdout));
    end

    [symbolsOut, xyzOut] = readXyzFile(xyzFile);
    if numel(symbolsOut) ~= N
        error('G_openbabel_optimize_geometry:atomCountChanged', ...
            ['Open Babel returned %d atoms, expected %d -- this should not happen with ' ...
             '''--noh'' passed; please report this.'], numel(symbolsOut), N);
    end
end

% -------------------------------------------------------------------------
function obabelPath = locateObabel()
% Checks PATH first (MATLAB's SYSTEM runs a non-interactive, non-login
% shell, so a conda/Homebrew install that only modifies an interactive
% shell's profile will usually NOT be found this way -- hence the
% explicit candidate list below, including the exact conda environment
% this function's own documentation recommends creating).
    [status, out] = system('which obabel');
    if status == 0 && ~isempty(strtrim(out))
        obabelPath = strtrim(out);
        return
    end
    home = getenv('HOME');
    candidates = { ...
        fullfile(home, 'miniforge3', 'envs', 'obabel_env', 'bin', 'obabel'), ...
        fullfile(home, 'miniforge3', 'bin', 'obabel'), ...
        fullfile(home, 'miniconda3', 'envs', 'obabel_env', 'bin', 'obabel'), ...
        fullfile(home, 'anaconda3', 'envs', 'obabel_env', 'bin', 'obabel'), ...
        '/usr/local/bin/obabel', ...
        '/opt/homebrew/bin/obabel'};
    for i = 1:numel(candidates)
        if isfile(candidates{i})
            obabelPath = candidates{i};
            return
        end
    end
    error('G_openbabel_optimize_geometry:obabelNotFound', ...
        ['Open Babel (''obabel'') was not found on PATH or in any common install location. ' ...
         'Install it, e.g.:\n  conda create -n obabel_env -c conda-forge openbabel -y\n' ...
         'then either add it to PATH or pass its full path via the ''ObabelPath'' option.']);
end

% -------------------------------------------------------------------------
function s = quoteForShell(pathStr)
% Wraps PATHSTR in single quotes for safe use inside a SYSTEM command
% (handles spaces; a path containing a literal single quote, extremely
% unlikely for a TEMPNAME-generated temp path, is not handled).
    s = ['''', pathStr, ''''];
end

% -------------------------------------------------------------------------
function writeMolfileV2000(filename, symbols, xyz, bondList)
% Writes an MDL Molfile (V2000): a header block, an atom block (element +
% Cartesian coordinates), and an explicit bond block (atom index pairs +
% bond order) -- the bond block is exactly why this format is used here
% instead of a bare .xyz, which carries no connectivity at all and would
% force Open Babel to re-perceive bonds from geometry alone.
    N = numel(symbols);
    if size(bondList, 2) < 3
        bondList = [bondList, ones(size(bondList, 1), 1)];
    end
    nBonds = size(bondList, 1);

    fid = fopen(filename, 'w');
    if fid < 0
        error('G_openbabel_optimize_geometry:cannotWrite', 'Cannot create temporary file: %s', filename);
    end
    fprintf(fid, 'Molecule\n');
    fprintf(fid, '  Generated by G_openbabel_optimize_geometry\n');
    fprintf(fid, '\n');
    fprintf(fid, '%3d%3d  0  0  0  0  0  0  0  0999 V2000\n', N, nBonds);
    for i = 1:N
        fprintf(fid, '%10.4f%10.4f%10.4f %-3s0  0  0  0  0  0  0  0  0  0  0  0\n', ...
            xyz(i,1), xyz(i,2), xyz(i,3), symbols{i});
    end
    for b = 1:nBonds
        order = bondList(b,3);
        if ~ismember(order, [1 2 3])
            order = 1;
        end
        fprintf(fid, '%3d%3d%3d  0  0  0  0\n', bondList(b,1), bondList(b,2), order);
    end
    fprintf(fid, 'M  END\n');
    fclose(fid);
end

% -------------------------------------------------------------------------
function [symbolsOut, xyzOut] = readXyzFile(filename)
% Reads a plain .xyz file (atom count, comment line, then one
% "Symbol X Y Z" line per atom) -- the format Open Babel was asked to
% write the minimized result in.
    fid = fopen(filename, 'r');
    if fid < 0
        error('G_openbabel_optimize_geometry:cannotRead', 'Cannot read file: %s', filename);
    end
    n = str2double(strtrim(fgetl(fid)));
    fgetl(fid);  % comment line, discarded
    symbolsOut = cell(n, 1);
    xyzOut = zeros(n, 3);
    for i = 1:n
        tok = strsplit(strtrim(fgetl(fid)));
        symbolsOut{i} = tok{1};
        xyzOut(i,:) = cellfun(@str2double, tok(2:4));
    end
    fclose(fid);
end
