function nm = G09_nmodes(filename, varargin)
% G09_NMODES  Extracts normal mode displacement vectors from a Gaussian 09 output file.
%
%   nm = G09_NMODES(filename)
%   nm = G09_NMODES(filename, 'modes', [5 10 20])
%
%   The format of the normal coordinate block is identical to G16.
%   G09 has only one Harmonic frequencies section per file.
%
%   OUTPUT  struct nm with fields:
%       .freq       [Nmodes x 1]          frequencies (cm-1)
%       .IR         [Nmodes x 1]          IR intensities (KM/Mole)
%       .Raman      [Nmodes x 1]          Raman activities (A^4/AMU), [] if absent
%       .redmass    [Nmodes x 1]          reduced masses (AMU)
%       .frcconst   [Nmodes x 1]          force constants (mDyne/A)
%       .symmetry   {Nmodes x 1 cell}     symmetry labels
%       .disp       [Natoms x 3 x Nmodes] displacement vectors
%       .Nmodes     int
%       .Natoms     int
%       .has_Raman  logical
%       .filename   char
%       .RamanFr        [Nmodes x K] per-incident-light-frequency Raman
%                       activities (A^4/AMU), [] if the file has no
%                       "RamAct Fr= N--" rows (i.e. no CPHF=RdFreq). Column
%                       1 is always the static (0 cm-1) value -- numerically
%                       identical to .Raman -- columns 2..K are the
%                       dynamic/pre-resonance values, one per incident
%                       wavelength actually requested, in the same order
%                       as .IncidentLight.
%       .IncidentLight  [1 x K]    incident light frequencies in cm-1, as
%                       printed on the "Incident light (cm**-1): ..." line
%                       ([] if absent). IncidentLight(1) is always 0.00.
%       .has_RamanFr    logical    true if .RamanFr is non-empty
%
%   Optional parameters also include:
%       'Lines'  - pre-read cell array of file lines (from G09_READ_LINES),
%                  to skip re-reading the file when it has already been
%                  read elsewhere (e.g. G09_READ_ALL). Default {} (read
%                  the file normally).
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename', @ischar);
addParameter(p, 'modes',    [],    @isnumeric);
addParameter(p, 'Lines',    {},    @iscell);
parse(p, filename, varargin{:});
mode_sel = p.Results.modes;

lines = p.Results.Lines;
if isempty(lines)
    lines = G09_read_lines(filename);
end
N     = numel(lines);

% Find section
sec_idx = find(~cellfun(@isempty, strfind(lines, 'and normal coordinates:')));
if isempty(sec_idx)
    error('G09_nmodes: "normal coordinates" section not found in %s', filename);
end
sec_start = sec_idx(end);

% -------------------------------------------------------------------------
% Parse — identical logic to G16_nmodes
% -------------------------------------------------------------------------
freqs    = [];
IRs      = [];
Ramans   = [];
redmass  = [];
frcconst = [];
symmetry = {};
disp_all = [];
Natoms_det = 0;
RamanFrAccum  = {};   % RamanFrAccum{n} = accumulated column for incident freq index n
IncidentLight = [];   % [1 x K] cm-1, parsed once from the section header

end_pat = '^\s*(-{20,}|Thermochemistry|Zero-point|Normal termination|Leave Link)';
k = sec_start + 1;

while k <= N
    ln = lines{k};
    if k > sec_start && ~isempty(regexp(ln, end_pat, 'once')), break; end

    if ~isempty(regexp(ln, '^\s*Frequencies\s*--', 'once'))
        if k > 1
            prev = strtrim(lines{k-1});
            if isempty(strfind(prev, '--')) && ~isempty(prev)
                for tok = strsplit(prev)
                    t = tok{1};
                    if ~isempty(t) && t(1) >= 'A' && t(1) <= 'z'
                        symmetry{end+1} = t; %#ok<AGROW>
                    end
                end
            end
        end
        freqs = [freqs; sscanf(parse_rhs(ln), '%f')]; %#ok<AGROW>
        k = k+1; continue
    end

    if ~isempty(regexp(ln, '^\s*Red\. masses\s*--', 'once'))
        redmass = [redmass; sscanf(parse_rhs(ln), '%f')]; %#ok<AGROW>
        k = k+1; continue
    end
    if ~isempty(regexp(ln, '^\s*Frc consts\s*--', 'once'))
        frcconst = [frcconst; sscanf(parse_rhs(ln), '%f')]; %#ok<AGROW>
        k = k+1; continue
    end
    if ~isempty(regexp(ln, '^\s*IR Inten\s*--', 'once'))
        IRs = [IRs; sscanf(parse_rhs(ln), '%f')]; %#ok<AGROW>
        k = k+1; continue
    end
    if ~isempty(regexp(ln, '^\s*Raman Activ\s*--', 'once'))
        Ramans = [Ramans; sscanf(parse_rhs(ln), '%f')]; %#ok<AGROW>
        k = k+1; continue
    end

    % ---------- Incident light (cm**-1): ... (CPHF=RdFreq header) ----------
    if isempty(IncidentLight) && ~isempty(regexp(ln, '^\s*Incident light\s*\(cm\*\*-1\):', 'once'))
        IncidentLight = sscanf(parse_rhs_colon(ln), '%f')';
        k = k+1; continue
    end

    % ---------- RamAct Fr= N-- (CPHF=RdFreq per-frequency Raman) ----------
    tok = regexp(ln, '^\s*RamAct Fr=\s*(\d+)\s*--', 'tokens', 'once');
    if ~isempty(tok)
        fridx = str2double(tok{1});
        vals  = sscanf(parse_rhs(ln), '%f');
        if numel(RamanFrAccum) < fridx
            RamanFrAccum{fridx} = []; %#ok<AGROW>
        end
        RamanFrAccum{fridx} = [RamanFrAccum{fridx}; vals]; %#ok<AGROW>
        k = k+1; continue
    end

    % Displacement block
    if ~isempty(regexp(ln, '^\s*Atom\s+AN\s+X', 'once'))
        n_so_far = numel(freqs);
        ncols    = n_so_far - size(disp_all, 3);
        if isempty(disp_all), ncols = n_so_far; end
        ncols = max(1, min(3, ncols));

        k = k+1;
        atom_disp = [];
        while k <= N
            ln2 = lines{k};
            tok2 = sscanf(ln2, '%f');
            if numel(tok2) < 2 + 3*ncols, break; end
            atom_disp(end+1, :) = tok2(3 : 2+3*ncols)'; %#ok<AGROW>
            k = k+1;
        end
        if isempty(atom_disp), continue; end

        Nat = size(atom_disp, 1);
        if Natoms_det == 0, Natoms_det = Nat; end

        for ci = 1:ncols
            col_xyz  = atom_disp(:, (ci-1)*3+1 : ci*3);
            disp_all = cat(3, disp_all, reshape(col_xyz, Nat, 3, 1));
        end
        continue
    end

    k = k+1;
end

Nmodes = numel(freqs);
if Nmodes == 0
    error('G09_nmodes: no modes read from %s', filename);
end

fix_vec = @(v) [v(1:min(end,Nmodes)); zeros(max(0,Nmodes-numel(v)),1)];
IRs      = fix_vec(IRs);
redmass  = fix_vec(redmass);
frcconst = fix_vec(frcconst);

has_Raman = numel(Ramans) == Nmodes;
if ~has_Raman, Ramans = []; end

% ---------- RamanFr: assemble the per-incident-frequency matrix ----------
if isempty(RamanFrAccum)
    RamanFr = [];
else
    K = numel(RamanFrAccum);
    RamanFr = NaN(Nmodes, K);
    for c = 1:K
        if ~isempty(RamanFrAccum{c})
            RamanFr(:, c) = fix_vec(RamanFrAccum{c});
        end
    end
end
has_RamanFr = ~isempty(RamanFr);

while numel(symmetry) < Nmodes, symmetry{end+1} = '?'; end
symmetry = symmetry(1:Nmodes);

% Optional mode selection
if ~isempty(mode_sel)
    idx      = mode_sel(mode_sel >= 1 & mode_sel <= Nmodes);
    freqs    = freqs(idx);
    IRs      = IRs(idx);
    if has_Raman, Ramans = Ramans(idx); end
    if has_RamanFr, RamanFr = RamanFr(idx, :); end
    redmass  = redmass(idx);
    frcconst = frcconst(idx);
    symmetry = symmetry(idx);
    if ~isempty(disp_all), disp_all = disp_all(:,:,idx); end
    Nmodes   = numel(freqs);
end

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
nm.RamanFr       = RamanFr;
nm.IncidentLight = IncidentLight;
nm.has_RamanFr   = has_RamanFr;

fprintf('G09_nmodes: %d modes, %d atoms — %s\n', Nmodes, Natoms_det, filename);

end  % G09_nmodes


function s = parse_rhs(ln)
idx = strfind(ln, '--');
s   = ln(idx(1)+2 : end);
end


% Right-hand side after the FIRST ':' (for the "Incident light (cm**-1):
% 0.00  12738.85" header line, whose numbers follow a colon, not a "--")
function s = parse_rhs_colon(ln)
idx = strfind(ln, ':');
if isempty(idx)
    s = ln;
else
    s = ln(idx(1)+1 : end);
end
end
