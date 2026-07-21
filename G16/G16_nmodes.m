function nm = G16_nmodes(filename, varargin)
% G16_NMODES  Extracts normal mode displacement vectors from a file
%             in Gaussian 16 .out/.log format.
%
%   nm = G16_NMODES(filename)
%   nm = G16_NMODES(filename, 'section', 'last')   % default
%   nm = G16_NMODES(filename, 'section', 'first')
%   nm = G16_NMODES(filename, 'modes', [5 10 15])  % load selected modes only
%
%   OUTPUT  struct nm with fields:
%       .freq        [Nmodes x 1]          frequencies (cm-1)
%       .IR          [Nmodes x 1]          IR intensities (KM/Mole)
%       .Raman       [Nmodes x 1]          Raman activities (A^4/AMU), [] if absent
%       .redmass     [Nmodes x 1]          reduced masses (AMU)
%       .frcconst    [Nmodes x 1]          force constants (mDyne/A)
%       .symmetry    {Nmodes x 1 cell}     symmetry labels
%       .disp        [Natoms x 3 x Nmodes] Cartesian displacement vectors
%       .Nmodes      int
%       .Natoms      int
%       .has_Raman   logical
%       .filename    char

% -------------------------------------------------------------------------
% Parse input arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename',  @ischar);
addParameter(p, 'section',   'last', @ischar);
addParameter(p, 'modes',     [],     @isnumeric);
addParameter(p, 'Lines',     {},     @iscell);
parse(p, filename, varargin{:});

sec_req  = lower(p.Results.section);
mode_sel = p.Results.modes;

% -------------------------------------------------------------------------
% Read file
% -------------------------------------------------------------------------
lines = p.Results.Lines;
if isempty(lines)
    if ~isfile(filename)
        error('G16_nmodes: file not found: %s', filename);
    end
    fid  = fopen(filename, 'r');
    raw  = fread(fid, '*char')';
    fclose(fid);
    lines = strsplit(raw, newline);
    G16_check_gaussian_match(lines, filename);
end
N = numel(lines);

% -------------------------------------------------------------------------
% Find the "and normal coordinates:" section
% -------------------------------------------------------------------------
idx_sec = find(~cellfun(@isempty, strfind(lines, 'and normal coordinates:')));
if isempty(idx_sec)
    error('G16_nmodes: "normal coordinates" section not found in %s', filename);
end

switch sec_req
    case 'last',  sec_start = idx_sec(end);
    case 'first', sec_start = idx_sec(1);
    otherwise
        error('G16_nmodes: section must be ''first'' or ''last''.');
end

% -------------------------------------------------------------------------
% Parse: linear scan
% -------------------------------------------------------------------------
% Logic: within the section, search cyclically for:
%   1) symmetry label row      (immediately before Frequencies --)
%   2) Frequencies --          vibrational frequencies
%   3) Red. masses --          reduced masses
%   4) Frc consts  --          force constants
%   5) IR Inten    --          IR intensities
%   6) [Raman Activ --]   optional
%   7) Atom AN X Y Z ...       read Natoms displacement rows
% Depolar, RamAct Fr=, Dep-P, Dep-U lines are skipped.

freqs    = [];
IRs      = [];
Ramans   = [];
redmass  = [];
frcconst = [];
symmetry = {};
disp_all = [];
Natoms_det = 0;

k = sec_start + 1;

while k <= N
    ln = lines{k};

    % ---------- end of section: exit loop ----------
    if ~isempty(regexp(ln, ...
            '^\s*(-{20,}|Thermochemistry|Zero-point|Normal termination|Leave Link)', ...
            'once'))
        break
    end

    % ---------- Frequencies -- ----------
    if ~isempty(regexp(ln, '^\s*Frequencies\s*--', 'once'))
        % The previous line (k-1) contains symmetry labels
        % es: "                      A                      A                      A"
        % Read labels from the previous line if it does not contain '--'
        if k > 1
            ln_prev = strtrim(lines{k-1});
            if isempty(strfind(ln_prev, '--')) && ~isempty(ln_prev)
                syms_raw = strsplit(ln_prev);
                syms_raw = syms_raw(~cellfun(@isempty, syms_raw));
                for ci = 1:numel(syms_raw)
                    % keep only tokens that begin with a letter (valid symmetry labels)
                    if ~isempty(regexp(syms_raw{ci}, '^[A-Za-z]', 'once'))
                        symmetry{end+1} = syms_raw{ci}; %#ok<AGROW>
                    end
                end
            end
        end
        vals = sscanf(parse_rhs(ln), '%f');
        freqs = [freqs; vals]; %#ok<AGROW>
        k = k+1; continue
    end

    % ---------- Red. masses -- ----------
    if ~isempty(regexp(ln, '^\s*Red\. masses\s*--', 'once'))
        vals = sscanf(parse_rhs(ln), '%f');
        redmass = [redmass; vals]; %#ok<AGROW>
        k = k+1; continue
    end

    % ---------- Frc consts -- ----------
    if ~isempty(regexp(ln, '^\s*Frc consts\s*--', 'once'))
        vals = sscanf(parse_rhs(ln), '%f');
        frcconst = [frcconst; vals]; %#ok<AGROW>
        k = k+1; continue
    end

    % ---------- IR Inten -- ----------
    if ~isempty(regexp(ln, '^\s*IR Inten\s*--', 'once'))
        vals = sscanf(parse_rhs(ln), '%f');
        IRs = [IRs; vals]; %#ok<AGROW>
        k = k+1; continue
    end

    % ---------- Raman Activ -- ----------
    if ~isempty(regexp(ln, '^\s*Raman Activ\s*--', 'once'))
        vals = sscanf(parse_rhs(ln), '%f');
        Ramans = [Ramans; vals]; %#ok<AGROW>
        k = k+1; continue
    end

    % ---------- Atom AN X Y Z ... → displacement block ----------
    if ~isempty(regexp(ln, '^\s*Atom\s+AN\s+X', 'once'))
        % Number of mode columns in this block:
        % = number of freqs read so far - size(disp_all,3)
        n_so_far = numel(freqs);
        if isempty(disp_all)
            ncols = n_so_far;
        else
            ncols = n_so_far - size(disp_all, 3);
        end
        if ncols < 1 || ncols > 3
            ncols = 3;  % fallback
        end

        k = k+1;
        atom_disp = [];

        while k <= N
            ln2 = lines{k};
            % Atom line: "  1  6  dx dy dz  [dx dy dz  [dx dy dz]]"
            % Minimo: iatom AN dx dy dz = 5 numeri
            tok2 = sscanf(ln2, '%f');
            if numel(tok2) < 2 + 3*ncols
                break
            end
            % tok2 = [iatom, AN, dx1,dy1,dz1, ...]
            row = tok2(3 : 2 + 3*ncols)';
            atom_disp(end+1, :) = row; %#ok<AGROW>
            k = k+1;
        end

        if isempty(atom_disp)
            continue
        end

        Nat = size(atom_disp, 1);
        if Natoms_det == 0
            Natoms_det = Nat;
        end

        % Append to disp_all [Nat x 3 x Nmodes]
        for ci = 1:ncols
            col_xyz = atom_disp(:, (ci-1)*3+1 : ci*3);
            disp_all = cat(3, disp_all, reshape(col_xyz, Nat, 3, 1));
        end
        continue
    end

    k = k+1;
end

% -------------------------------------------------------------------------
% Check and align vectors
% -------------------------------------------------------------------------
Nmodes = numel(freqs);
if Nmodes == 0
    error('G16_nmodes: no modes read from %s', filename);
end

fix_vec = @(v) [v(1:min(end,Nmodes)); zeros(max(0,Nmodes-numel(v)), 1)];
IRs      = fix_vec(IRs);
redmass  = fix_vec(redmass);
frcconst = fix_vec(frcconst);

has_Raman = numel(Ramans) == Nmodes;
if ~has_Raman
    Ramans = [];
end

if numel(symmetry) < Nmodes
    symmetry(end+1:Nmodes) = {'?'};
end
symmetry = symmetry(1:Nmodes);

% -------------------------------------------------------------------------
% Selezione modi optional
% -------------------------------------------------------------------------
if ~isempty(mode_sel)
    mode_sel = mode_sel(mode_sel >= 1 & mode_sel <= Nmodes);
    freqs    = freqs(mode_sel);
    IRs      = IRs(mode_sel);
    if has_Raman, Ramans = Ramans(mode_sel); end
    redmass  = redmass(mode_sel);
    frcconst = frcconst(mode_sel);
    symmetry = symmetry(mode_sel);
    if ~isempty(disp_all)
        disp_all = disp_all(:, :, mode_sel);
    end
    Nmodes = numel(freqs);
end

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
nm.freq      = freqs;
nm.IR        = IRs;
nm.Raman     = Ramans;
nm.redmass   = redmass;
nm.frcconst  = frcconst;
nm.symmetry  = symmetry;
nm.disp      = disp_all;
nm.Nmodes    = Nmodes;
nm.Natoms    = Natoms_det;
nm.has_Raman = has_Raman;
nm.filename  = filename;

fprintf('G16_nmodes: %d modes, %d atoms read from %s\n', Nmodes, Natoms_det, filename);

end  % G16_nmodes


% =========================================================================
%  Local function: right-hand side after '--'
% =========================================================================
function s = parse_rhs(ln)
idx = strfind(ln, '--');
if isempty(idx)
    s = ln;
else
    s = ln(idx(1)+2 : end);
end
end
