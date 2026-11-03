function S = G16_scan_read(filename, varargin)
% G16_SCAN_READ  Reads a relaxed PES scan (Opt=ModRedundant with an 'S' coordinate) point by point
%
%   S = G16_SCAN_READ(FILENAME) parses a Gaussian 16 relaxed scan and
%   returns, for every scan point, the converged constrained geometry, its
%   SCF energy, the scanned coordinate(s) measured on that geometry and the
%   energy relative to the lowest point. The scan job is located
%   automatically when the file holds several jobs chained with --Link1--.
%
%   S = G16_SCAN_READ(FILENAME, 'Name', Value, ...)
%
%   NAME-VALUE OPTIONS:
%     'job'   - 'auto' (default) | integer: which job of a --Link1-- chain
%               to read (1-based); 'auto' takes the first relaxed scan
%     'plot'  - false (default) | true: plot dE (kcal/mol) versus the
%               scanned coordinate (1-D scan) or the point index
%     'Lines' - pre-read file lines (see G16_read_all)
%
%   OUTPUT:
%     S - struct with fields:
%       .coords        struct array of scanned coordinates: .type ('B','A',
%                      'D'), .atoms (1-based), .label (e.g. 'D(4,1,3,9)'),
%                      .nsteps, .stepsize, .units
%       .frozen        cell array of frozen ('F') coordinates, as labels
%       .value         [Npts x Ncoord] scanned coordinate(s) measured on
%                      each geometry (dihedrals unwrapped along the scan)
%       .E             [Npts x 1] SCF energy (Hartree)
%       .dE_kcal, .dE_kJ  [Npts x 1] energy relative to the lowest point
%       .i_min         index of the lowest point
%       .point         [Npts x 1] scan-point number
%       .converged     [Npts x 1] logical ('Stationary point found')
%       .opt_steps     [Npts x 1] optimisation steps used per point
%       .n_points_expected  number of points announced by Gaussian
%       .mol           struct array, one per point, compatible with
%                      G16_get_bond_length / G16_get_angle /
%                      G16_get_dihedral / G16_draw_molecule
%       .xyz           [Natoms x 3 x Npts] coordinates (Angstrom)
%       .symbols, .Z, .Natoms, .orientation, .method
%       .completed     true if the job ended with Normal termination
%       .status        termination message
%       .job           job index used
%       .filename
%
%   NOTE: Gaussian prints the dipole moment and atomic charges of a
%   relaxed scan only for the first and last geometry, so they are not
%   returned per point; run single points on S.mol(k) if needed.
%
%   EXAMPLE:
%     S = G16_scan_read('NMA_scanCN.log', 'plot', true);
%     [~, phi] = G16_get_dihedral(S.mol(13), [4 1 3 9]);
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

p = inputParser;
addRequired(p,  'filename', @ischar);
addParameter(p, 'job',   'auto', @(x) ischar(x) || isnumeric(x));
addParameter(p, 'plot',  false,  @islogical);
addParameter(p, 'Lines', {},     @iscell);
parse(p, filename, varargin{:});
lines = p.Results.Lines;
if isempty(lines)
    if exist(filename, 'file') ~= 2
        error('G16_scan_read: file not found: %s', filename);
    end
    fid = fopen(filename, 'r');
    raw = fread(fid, '*char')';
    fclose(fid);
    raw = strrep(raw, sprintf('\r'), '');
    lines = strsplit(raw, sprintf('\n'));
    G16_check_gaussian_match(lines, filename);
end
HA2KCAL = 627.5095;
HA2KJ   = 2625.4996;

% ---------- job selection --------------------------------------------
jobs = split_jobs(lines);
if ischar(p.Results.job)
    jidx = 0;
    for n = 1:size(jobs, 1)
        if any(~cellfun(@isempty, regexp(lines(jobs(n,1):jobs(n,2)), 'on scan point', 'once')))
            jidx = n; break
        end
    end
    if jidx == 0
        error('G16_scan_read: no relaxed scan (''on scan point'') found in %s', filename);
    end
else
    jidx = round(p.Results.job);
    if jidx < 1 || jidx > size(jobs, 1)
        error('G16_scan_read: job %d out of range [1, %d].', jidx, size(jobs, 1));
    end
end
a = jobs(jidx, 1);  b = jobs(jidx, 2);

% ---------- scanned / frozen coordinates (ModRedundant echo) ---------
coords = struct('type', {}, 'atoms', {}, 'label', {}, 'nsteps', {}, ...
                'stepsize', {}, 'units', {});
frozen = {};
for k = a:b
    if isempty(strfind(lines{k}, 'ModRedundant input section has been read')), continue; end
    for j = k+1:b
        tok = regexp(lines{j}, ...
            '^\s*([BADLbadl])\s+((?:\d+\s+){2,4})([SFsf])(?=\s|$)(\s+\d+\s+[-\d.]+)?', ...
            'tokens', 'once');
        if isempty(tok), break; end
        % The trailing "nsteps stepsize" span is captured as a single group
        % (parsed below with sscanf): a capturing group nested inside an
        % optional non-capturing group, e.g. (?:\s+(\d+)\s+([-\d.]+))?, is
        % silently dropped from MATLAB's 'tokens' output even when it does
        % participate (unlike Octave, where this toolbox was first tried).
        tok(end+1:4) = {''};
        ctype = upper(tok{1});
        atoms = sscanf(tok{2}, '%d')';
        label = sprintf('%s(%s)', ctype, strjoin(arrayfun(@num2str, atoms, ...
                        'UniformOutput', false), ','));
        if upper(tok{3}) == 'S'
            c.type = ctype; c.atoms = atoms; c.label = label;
            sv = sscanf(tok{4}, '%d %f')';
            if numel(sv) == 2
                c.nsteps = sv(1); c.stepsize = sv(2);
            else
                c.nsteps = NaN; c.stepsize = NaN;
            end
            if ctype == 'B', c.units = 'Angstrom'; else, c.units = 'deg'; end
            coords(end+1) = c; %#ok<AGROW>
        else
            frozen{end+1} = label; %#ok<AGROW>
        end
    end
    break
end

% ---------- markers ---------------------------------------------------
ori_label = 'Input orientation';
for k = a:b
    if ~isempty(strfind(lines{k}, 'Standard orientation')), ori_label = 'Standard orientation'; break; end
end
scf_idx = []; scf_val = []; method = '';
ori_idx = []; mark_idx = []; mark_pt = []; stat_idx = [];
steps = containers.Map('KeyType', 'double', 'ValueType', 'double');
n_expected = 0;
for k = a:b
    ln = lines{k};
    tok = regexp(ln, 'SCF Done:\s+E\((\S+)\)\s*=\s*([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok)
        scf_idx(end+1) = k; scf_val(end+1) = str2double(tok{2}); method = tok{1}; %#ok<AGROW>
        continue
    end
    if ~isempty(strfind(ln, ori_label)) && ~isempty(regexp(ln, 'orientation\s*:', 'once'))
        ori_idx(end+1) = k; %#ok<AGROW>
        continue
    end
    tok = regexp(ln, 'on scan point\s+(\d+)\s+out of\s+(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        pt = str2double(tok{1});
        mark_idx(end+1) = k; mark_pt(end+1) = pt; %#ok<AGROW>
        n_expected = str2double(tok{2});
        ts = regexp(ln, 'Step number\s+(\d+)', 'tokens', 'once');
        if ~isempty(ts)
            ns = str2double(ts{1});
            if ~isKey(steps, pt) || steps(pt) < ns, steps(pt) = ns; end
        end
        continue
    end
    if ~isempty(strfind(ln, 'Stationary point found'))
        stat_idx(end+1) = k; %#ok<AGROW>
    end
end
if isempty(mark_idx)
    error('G16_scan_read: job %d of %s is not a relaxed scan.', jidx, filename);
end

% point ends: (line, point, converged)
ends = zeros(0, 3);
for k = stat_idx
    j = last_before(mark_idx, k);
    if isempty(j), continue; end
    pt = mark_pt(j);
    if isempty(ends) || ~any(ends(:,2) == pt)
        ends(end+1, :) = [k, pt, 1]; %#ok<AGROW>
    end
end
if isempty(ends) || ~any(ends(:,2) == mark_pt(end))
    ends(end+1, :) = [b + 1, mark_pt(end), 0];   % job stopped inside this point
end

E = []; pts = []; conv = []; XYZc = {}; Zref = [];
for r = 1:size(ends, 1)
    js = last_before(scf_idx, ends(r,1));
    jo = last_before(ori_idx, ends(r,1));
    if isempty(js) || isempty(jo), continue; end
    [Z, xyz] = read_orientation_block(lines, ori_idx(jo));
    if isempty(Zref), Zref = Z; end
    E(end+1, 1)    = scf_val(js);   %#ok<AGROW>
    pts(end+1, 1)  = ends(r,2);     %#ok<AGROW>
    conv(end+1, 1) = ends(r,3);     %#ok<AGROW>
    XYZc{end+1} = xyz;              %#ok<AGROW>
end
if isempty(E)
    error('G16_scan_read: no scan point could be read from %s', filename);
end
XYZ = cat(3, XYZc{:});
Npts = numel(E);
symbols = arrayfun(@z2sym, Zref, 'UniformOutput', false);

value = nan(Npts, numel(coords));
for c = 1:numel(coords)
    for n = 1:Npts
        value(n, c) = measure(XYZ(:, :, n), coords(c).type, coords(c).atoms);
    end
    if any(coords(c).type == 'DL')
        value(:, c) = rad2deg_(unwrap(deg2rad_(value(:, c))));
    end
end

[~, i_min] = min(E);
dE = E - E(i_min);
[completed, status] = job_status(lines, a, b);
if any(~conv)
    warning('G16_scan_read:notConverged', '%d scan point(s) did not converge in %s (%s).', ...
            sum(~conv), filename, status);
end

mol = struct('symbols', {}, 'xyz', {}, 'Z', {}, 'Natoms', {}, 'step', {}, ...
             'n_steps', {}, 'orientation', {}, 'filename', {});
for n = 1:Npts
    mol(n).symbols = symbols(:); mol(n).xyz = XYZ(:, :, n); mol(n).Z = Zref(:);
    mol(n).Natoms = numel(Zref); mol(n).step = n; mol(n).n_steps = Npts;
    mol(n).orientation = ori_label; mol(n).filename = filename;
end
opt_steps = zeros(Npts, 1);
for n = 1:Npts
    if isKey(steps, pts(n)), opt_steps(n) = steps(pts(n)); end
end

S.coords  = coords;
S.frozen  = frozen;
S.value   = value;
S.E       = E;
S.dE_kcal = dE * HA2KCAL;
S.dE_kJ   = dE * HA2KJ;
S.i_min   = i_min;
S.point   = pts;
S.converged = logical(conv);
S.opt_steps = opt_steps;
S.n_points_expected = n_expected;
S.mol     = mol;
S.xyz     = XYZ;
S.symbols = symbols(:);
S.Z       = Zref(:);
S.Natoms  = numel(Zref);
S.orientation = ori_label;
S.method  = method;
S.completed = completed;
S.status  = status;
S.job     = jidx;
S.filename = filename;

if isempty(coords), lab = '(no S coordinate echoed)';
else, lab = strjoin({coords.label}, ', '); end
fprintf('\n── G16_scan_read: %s (job %d) ──\n', filename, jidx);
fprintf('  Scanned: %s   points read: %d/%d   %s\n', lab, Npts, n_expected, status);
fprintf('  E(min) = %.8f Ha at point %d;  dE(max) = %.3f kcal/mol\n\n', ...
        E(i_min), pts(i_min), max(S.dE_kcal));

if p.Results.plot
    figure('Color', 'w');
    if numel(coords) == 1
        x = value(:, 1); xl = sprintf('%s (%s)', coords(1).label, coords(1).units);
    else
        x = pts; xl = 'Scan point';
    end
    plot(x, S.dE_kcal, 'o-', 'MarkerSize', 4); hold on
    if any(~S.converged)
        plot(x(~S.converged), S.dE_kcal(~S.converged), 'rx', 'MarkerSize', 9);
    end
    xlabel(xl); ylabel('\DeltaE (kcal/mol)');
    title(sprintf('Relaxed scan - %s', method)); grid on; box on
end
end % G16_scan_read


% =========================================================================
function jobs = split_jobs(lines)
% One row [first last] per job; --Link1-- jobs start with 'Entering Link 1'
% (the Freq step Gaussian appends to Opt+Freq does not, and stays with its Opt).
st = find(~cellfun(@isempty, strfind(lines, 'Entering Link 1')));
if isempty(st), jobs = [1, numel(lines)]; return; end
st(1) = 1;
jobs = [st(:), [st(2:end)' - 1; numel(lines)]];
end

function [ok, msg] = job_status(lines, a, b)
ok = false; msg = 'No termination line (job still running or file truncated)';
for k = b:-1:a
    if ~isempty(strfind(lines{k}, 'Normal termination'))
        ok = true; msg = 'Normal termination'; return
    end
    if ~isempty(strfind(lines{k}, 'Error termination'))
        msg = 'Error termination';
        for j = k-1:-1:max(a, k-6)
            s = strtrim(lines{j});
            if ~isempty(s) && isempty(regexp(s, '^(Error|Job cpu|Elapsed|File lengths)', 'once'))
                msg = ['Error termination: ' s]; break
            end
        end
        return
    end
end
end

function j = last_before(idx, k)
% index (into idx) of the largest element < k, [] if none
j = find(idx < k, 1, 'last');
end

function [Z, xyz] = read_orientation_block(lines, k)
Z = []; xyz = zeros(0, 3);
for j = k+5:numel(lines)
    s = strtrim(lines{j});
    if ~isempty(s) && all(s == '-'), break; end
    v = sscanf(s, '%f')';
    if numel(v) == 6
        Z(end+1, 1) = v(2); xyz(end+1, :) = v(4:6); %#ok<AGROW>
    end
end
end

function v = measure(xyz, ctype, at)
switch ctype
    case 'B'
        v = norm(xyz(at(1),:) - xyz(at(2),:));
    case 'A'
        u = xyz(at(1),:) - xyz(at(2),:); w = xyz(at(3),:) - xyz(at(2),:);
        v = acos(max(-1, min(1, dot(u, w) / (norm(u) * norm(w))))) * 180 / pi;
    case {'D', 'L'}
        b0 = xyz(at(1),:) - xyz(at(2),:); b1 = xyz(at(3),:) - xyz(at(2),:);
        b2 = xyz(at(4),:) - xyz(at(3),:);
        b1 = b1 / norm(b1);
        u = b0 - dot(b0, b1) * b1; w = b2 - dot(b2, b1) * b1;
        v = atan2(dot(cross(b1, u), w), dot(u, w)) * 180 / pi;
    otherwise
        v = NaN;
end
end

function r = deg2rad_(d), r = d * pi / 180; end
function d = rad2deg_(r), d = r * 180 / pi; end

function s = z2sym(z)
persistent T
if isempty(T)
    T = {'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S', ...
         'Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga', ...
         'Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd', ...
         'Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm', ...
         'Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os', ...
         'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa', ...
         'U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg', ...
         'Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'};
end
s = T{z};
end
