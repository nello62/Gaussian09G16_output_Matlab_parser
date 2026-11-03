function I = G16_irc_read(filename, varargin)
% G16_IRC_READ  Reads an IRC calculation point by point, ordered reverse -> TS -> forward
%
%   I = G16_IRC_READ(FILENAME) parses a Gaussian 16 IRC and returns every
%   point Gaussian completed: geometry, energy, signed reaction coordinate
%   (negative on the REVERSE branch, 0 at the TS) and energy relative to
%   the TS. Works on a stand-alone IRC job, on an IRC that is one job of a
%   --Link1-- chain (e.g. TS optimisation + Freq + IRC), and on an IRC
%   that stopped with an error.
%
%   I = G16_IRC_READ(FILENAME, 'Name', Value, ...)
%
%   NAME-VALUE OPTIONS:
%     'job'   - 'auto' (default) | integer: which job of a --Link1-- chain
%               to read (1-based); 'auto' takes the first IRC
%     'plot'  - false (default) | true: plot dE (kcal/mol, relative to the
%               TS) versus the reaction coordinate
%     'Lines' - pre-read file lines (see G16_read_all)
%
%   OUTPUT:
%     I - struct with fields:
%       .rxcoord      [Npts x 1] reaction coordinate (sqrt(amu)*bohr)
%       .E            [Npts x 1] energy (Hartree), NaN if not printed
%       .dE_kcal, .dE_kJ  [Npts x 1] energy relative to the TS
%       .E_source     cell: 'chk' (TS energy read from checkpoint), 'scf'
%                     (SCF at the point) or 'hpc_corrected' (HPC corrector
%                     end-point energy, DWI-interpolated)
%       .direction    cell: 'REVERSE' | 'TS' | 'FORWARD'
%       .point        [Npts x 1] Gaussian point number on its branch
%       .i_ts         index of the TS
%       .n_forward, .n_reverse   points completed on each branch
%       .mol          struct array, one per point, compatible with
%                     G16_get_bond_length / G16_get_angle /
%                     G16_get_dihedral / G16_draw_molecule
%       .xyz          [Natoms x 3 x Npts] coordinates (Angstrom)
%       .symbols, .Z, .Natoms
%       .max_points   MaxPoints per path requested
%       .step_size    e.g. '0.100 bohr'
%       .scheme       integration scheme (HPC, LQA, ...)
%       .completed    true if the job ended with Normal termination
%       .status       termination message
%       .job          job index used
%       .filename
%
%   EXAMPLE:
%     I = G16_irc_read('NMA_TS_IRC_CN.log', 'plot', true);
%     [~, phi] = G16_get_dihedral(I.mol(1), [4 1 3 9]);   % reverse end
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
        error('G16_irc_read: file not found: %s', filename);
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
        if any(~cellfun(@isempty, strfind(lines(jobs(n,1):jobs(n,2)), 'IRC-IRC-IRC')))
            jidx = n; break
        end
    end
    if jidx == 0
        error('G16_irc_read: no IRC calculation found in %s', filename);
    end
else
    jidx = round(p.Results.job);
    if jidx < 1 || jidx > size(jobs, 1)
        error('G16_irc_read: job %d out of range [1, %d].', jidx, size(jobs, 1));
    end
end
a = jobs(jidx, 1);  b = jobs(jidx, 2);

% ---------- markers ---------------------------------------------------
scf_idx = []; scf_val = []; ori_idx = [];
corr_idx = []; corr_val = [];
dir_idx = []; dir_pt = []; dir_name = {};
pt_idx = []; pt_num = []; pt_path = [];
max_points = NaN; step_size = ''; scheme = ''; E_chk = NaN;
for k = a:b
    ln = lines{k};
    tok = regexp(ln, 'SCF Done:\s+E\((\S+)\)\s*=\s*([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), scf_idx(end+1) = k; scf_val(end+1) = str2double(tok{2}); continue; end %#ok<AGROW>
    if ~isempty(regexp(ln, '(Standard|Input) orientation\s*:', 'once')), ori_idx(end+1) = k; continue; end %#ok<AGROW>
    tok = regexp(ln, 'Corrected End Point Energy\s*=\s*([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok), corr_idx(end+1) = k; corr_val(end+1) = str2double(tok{1}); continue; end %#ok<AGROW>
    tok = regexp(ln, 'Point Number:\s+(\d+)\s+Path Number:\s+(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        pt_idx(end+1) = k; pt_num(end+1) = str2double(tok{1}); pt_path(end+1) = str2double(tok{2}); %#ok<AGROW>
        continue
    end
    tok = regexp(ln, 'Point Number\s+(\d+)\s+in\s+(FORWARD|REVERSE)\s+path direction', 'tokens', 'once');
    if ~isempty(tok)
        dir_idx(end+1) = k; dir_pt(end+1) = str2double(tok{1}); dir_name{end+1} = upper(tok{2}); %#ok<AGROW>
        continue
    end
    tok = regexp(ln, 'Energy From Chk\s*=\s*([-\d.]+)', 'tokens', 'once');
    if ~isempty(tok) && isnan(E_chk), E_chk = str2double(tok{1}); continue; end
    tok = regexp(ln, 'Maximum points per path\s*=\s*(\d+)', 'tokens', 'once');
    if ~isempty(tok) && isnan(max_points), max_points = str2double(tok{1}); continue; end
    tok = regexp(ln, 'Step size\s*=\s*([\d.]+)\s*(\S+)', 'tokens', 'once');
    if ~isempty(tok) && isempty(step_size), step_size = [tok{1} ' ' tok{2}]; continue; end
    tok = regexp(ln, 'Integration scheme\s*=\s*(\S+)', 'tokens', 'once');
    if ~isempty(tok) && isempty(scheme), scheme = tok{1}; end
end
if isempty(pt_idx)
    error('G16_irc_read: no IRC point found in job %d of %s', jidx, filename);
end

% ---------- one row per point ----------------------------------------
R = struct('direction', {}, 'point', {}, 'rx', {}, 'E', {}, 'src', {}, 'Z', {}, 'xyz', {});
seen = {};
prev = a;
for q = 1:numel(pt_idx)
    k = pt_idx(q);  pn = pt_num(q);
    if pn == 0
        if any(strcmp({R.direction}, 'TS')), prev = k; continue; end
        jo = last_before(ori_idx, k);
        if ~isempty(jo)
            [Z, xyz] = read_orientation_block(lines, ori_idx(jo));
        else
            [Z, xyz] = read_checkpoint_structure(lines, a, k);
        end
        if ~isnan(E_chk)
            E = E_chk; src = 'chk';
        else
            js = last_before(scf_idx, k);
            if isempty(js), E = NaN; src = 'none'; else, E = scf_val(js); src = 'scf'; end
        end
        R(end+1) = struct('direction', 'TS', 'point', 0, 'rx', 0, 'E', E, 'src', src, ...
                          'Z', Z, 'xyz', xyz); %#ok<AGROW>
        prev = k;
        continue
    end
    jd = last_before(dir_idx, k);
    if ~isempty(jd)
        direction = dir_name{jd};
    elseif pt_path(q) == 1
        direction = 'FORWARD';
    else
        direction = 'REVERSE';
    end
    key = sprintf('%s%d', direction, pn);
    if any(strcmp(seen, key)), prev = k; continue; end
    seen{end+1} = key; %#ok<AGROW>

    [Z, xyz] = read_cartesian_block(lines, k + 1);
    if isempty(Z)
        jo = last_before(ori_idx, k);
        if ~isempty(jo), [Z, xyz] = read_orientation_block(lines, ori_idx(jo)); end
    end
    rx = NaN;
    for j = k+1:min(b, k+40)
        tok = regexp(lines{j}, 'NET REACTION COORDINATE UP TO THIS POINT\s*=\s*([-\d.]+)', 'tokens', 'once');
        if ~isempty(tok), rx = str2double(tok{1}); break; end
    end
    jc = last_before(corr_idx, k);
    js = last_before(scf_idx, k);
    if ~isempty(jc) && corr_idx(jc) > prev
        E = corr_val(jc); src = 'hpc_corrected';
    elseif ~isempty(js) && scf_idx(js) > prev
        E = scf_val(js); src = 'scf';
    else
        E = NaN; src = 'none';
    end
    sgn = 1; if strcmp(direction, 'REVERSE'), sgn = -1; end
    if ~isempty(Z)
        R(end+1) = struct('direction', direction, 'point', pn, 'rx', sgn * abs(rx), ...
                          'E', E, 'src', src, 'Z', Z, 'xyz', xyz); %#ok<AGROW>
    end
    prev = k;
end

% ---------- order: reverse (far -> near), TS, forward (near -> far) ----
dirs = {R.direction};  pnum = [R.point];
ir = find(strcmp(dirs, 'REVERSE'));  [~, o] = sort(pnum(ir), 'descend');  ir = ir(o);
it = find(strcmp(dirs, 'TS'));
iF = find(strcmp(dirs, 'FORWARD'));  [~, o] = sort(pnum(iF), 'ascend');   iF = iF(o);
R = R([ir, it, iF]);
Npts = numel(R);
E = [R.E]';
if isempty(it)
    warning('G16_irc_read:noTS', ['TS (point 0) geometry not found in %s; energies are ' ...
            'given relative to the highest point instead.'], filename);
    [~, i_ts] = max(E);
else
    i_ts = numel(ir) + 1;
end
dE = E - E(i_ts);
Zref = R(1).Z;
symbols = arrayfun(@z2sym, Zref, 'UniformOutput', false);
XYZ = cat(3, R.xyz);
[completed, status] = job_status(lines, a, b);

mol = struct('symbols', {}, 'xyz', {}, 'Z', {}, 'Natoms', {}, 'step', {}, ...
             'n_steps', {}, 'orientation', {}, 'filename', {});
for n = 1:Npts
    mol(n).symbols = symbols(:); mol(n).xyz = R(n).xyz; mol(n).Z = Zref(:);
    mol(n).Natoms = numel(Zref); mol(n).step = n; mol(n).n_steps = Npts;
    mol(n).orientation = 'IRC path'; mol(n).filename = filename;
end

I.rxcoord   = [R.rx]';
I.E         = E;
I.dE_kcal   = dE * HA2KCAL;
I.dE_kJ     = dE * HA2KJ;
I.E_source  = {R.src}';
I.direction = {R.direction}';
I.point     = [R.point]';
I.i_ts      = i_ts;
I.n_forward = numel(iF);
I.n_reverse = numel(ir);
I.mol       = mol;
I.xyz       = XYZ;
I.symbols   = symbols(:);
I.Z         = Zref(:);
I.Natoms    = numel(Zref);
I.max_points = max_points;
I.step_size = step_size;
I.scheme    = scheme;
I.completed = completed;
I.status    = status;
I.job       = jidx;
I.filename  = filename;

if isempty(it), tslab = 'no TS'; else, tslab = 'TS'; end
fprintf('\n── G16_irc_read: %s (job %d) ──\n', filename, jidx);
fprintf('  Points: %d reverse + %s + %d forward   (MaxPoints=%d, step %s, %s)\n', ...
        I.n_reverse, tslab, I.n_forward, max_points, step_size, scheme);
fprintf('  %s\n', status);
if ~completed
    fprintf('  Note: the IRC stopped before the end; all completed points are returned.\n');
end
fprintf('\n');

if p.Results.plot
    figure('Color', 'w');
    plot(I.rxcoord, I.dE_kcal, 'o-', 'MarkerSize', 3); hold on
    plot(I.rxcoord(i_ts), I.dE_kcal(i_ts), 'p', 'MarkerSize', 12, ...
         'MarkerFaceColor', [0.85 0.1 0.2], 'MarkerEdgeColor', [0.85 0.1 0.2]);
    xlabel('Reaction coordinate (amu^{1/2} bohr)   \leftarrow reverse | forward \rightarrow');
    ylabel('\DeltaE relative to TS (kcal/mol)');
    if completed, title('IRC'); else, title('IRC (incomplete)'); end
    grid on; box on
end
end % G16_irc_read


% =========================================================================
function jobs = split_jobs(lines)
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

function [Z, xyz] = read_cartesian_block(lines, k)
% 5-column 'Center / Atomic Number / X Y Z' rows (IRC 'CURRENT STRUCTURE')
Z = []; xyz = zeros(0, 3); started = false;
for j = k:numel(lines)
    v = sscanf(strtrim(lines{j}), '%f')';
    if numel(v) == 5 && v(1) == round(v(1)) && v(2) == round(v(2))
        started = true; Z(end+1, 1) = v(2); xyz(end+1, :) = v(3:5); %#ok<AGROW>
        continue
    end
    if started || j - k > 8, break; end
end
end

function [Z, xyz] = read_checkpoint_structure(lines, a, b)
% 'Structure from the checkpoint file' echo: lines 'C,0,x,y,z' or 'C,x,y,z'
Z = []; xyz = zeros(0, 3);
for k = a:b
    if isempty(strfind(lines{k}, 'Structure from the checkpoint file')), continue; end
    for j = k+1:b
        parts = strtrim(strsplit(strtrim(lines{j}), ','));
        ok = (numel(parts) == 4 || numel(parts) == 5) && ~isempty(parts{1}) && isletter(parts{1}(1));
        if ok
            v = str2double(parts(end-2:end));
            zz = sym2z(regexprep(parts{1}, '[^A-Za-z]', ''));
            ok = all(~isnan(v)) && ~isnan(zz);
        end
        if ok
            Z(end+1, 1) = zz; xyz(end+1, :) = v; %#ok<AGROW>
        elseif ~isempty(Z)
            break
        end
    end
    if ~isempty(Z), return; end
end
end

function s = z2sym(z)
T = sym_table();
s = T{z};
end

function z = sym2z(s)
T = sym_table();
z = find(strcmpi(T, s), 1);
if isempty(z), z = NaN; end
end

function T = sym_table()
persistent TT
if isempty(TT)
    TT = {'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S', ...
          'Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga', ...
          'Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd', ...
          'Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm', ...
          'Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os', ...
          'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa', ...
          'U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg', ...
          'Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'};
end
T = TT;
end
