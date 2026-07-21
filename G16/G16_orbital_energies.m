function oe = G16_orbital_energies(filename, varargin)
% G16_ORBITAL_ENERGIES  Extracts molecular orbital energies (HOMO/LUMO and
%                        the full occupied/virtual spectrum) from a
%                        Gaussian 09 output file.
%
%   oe = G16_ORBITAL_ENERGIES(filename)
%   oe = G16_ORBITAL_ENERGIES(filename, 'step', 'last')
%   oe = G16_ORBITAL_ENERGIES(filename, 'step', N)
%
%   Reads the "Alpha  occ. eigenvalues --" / "Alpha virt. eigenvalues --"
%   (and, for open-shell calculations, the matching "Beta" lines) printed
%   by Gaussian's population analysis. These blocks repeat once per SCF
%   calculation in the file (e.g. once per optimisation step); 'step'
%   selects which block to report on, exactly like G16_ENERGY/G16_STRUCTURE.
%
%   Optional parameters (Name-Value):
%       'step'   - 'last' (default) | 'first' | integer index
%
%   OUTPUT  struct oe with fields (all energies in Hartree unless noted):
%       .alpha_occ      [Nocc x 1]   occupied alpha orbital energies
%       .alpha_virt     [Nvirt x 1]  virtual alpha orbital energies
%       .beta_occ       [Nocc x 1]   occupied beta orbital energies ([] if closed-shell)
%       .beta_virt      [Nvirt x 1]  virtual beta orbital energies ([] if closed-shell)
%       .has_beta       logical      true for open-shell (UHF/UKS) calculations
%       .HOMO           double       highest occupied orbital energy
%       .LUMO           double       lowest unoccupied orbital energy
%       .gap            double       LUMO - HOMO
%       .HOMO_alpha     double       alpha HOMO (== .HOMO for closed-shell)
%       .LUMO_alpha     double       alpha LUMO (== .LUMO for closed-shell)
%       .HOMO_beta      double       beta HOMO ([] if closed-shell)
%       .LUMO_beta      double       beta LUMO ([] if closed-shell)
%       .HOMO_eV .LUMO_eV .gap_eV    same three quantities, in eV
%       .step           int          block index actually used
%       .Nsteps         int          number of eigenvalue blocks found in the file
%       .filename       char
%
%   Example:
%       oe = G16_orbital_energies('V_E00t.out');
%       fprintf('HOMO-LUMO gap = %.3f eV\n', oe.gap_eV);

% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename', @ischar);
addParameter(p, 'step',     'last', @(x) ischar(x) || isnumeric(x));
parse(p, filename, varargin{:});
step_req = p.Results.step;

if ~isfile(filename)
    error('G16_orbital_energies: file not found: %s', filename);
end
fid  = fopen(filename, 'r');
raw  = fread(fid, '*char')';
fclose(fid);
lines = strsplit(raw, newline);
G16_check_gaussian_match(lines, filename);
N     = numel(lines);

% -------------------------------------------------------------------------
% Scan the file, grouping consecutive eigenvalue lines into blocks. A new
% block starts whenever an "Alpha  occ." line is seen right after a
% "virt." line (Alpha or Beta) — i.e. a fresh population analysis.
% -------------------------------------------------------------------------
blocks = struct('alpha_occ', {}, 'alpha_virt', {}, 'beta_occ', {}, 'beta_virt', {});
cur_alpha_occ = []; cur_alpha_virt = []; cur_beta_occ = []; cur_beta_virt = [];
last_kind = '';

for k = 1:N
    ln = lines{k};

    if contains(ln, 'Alpha  occ. eigenvalues') || contains(ln, 'Alpha occ. eigenvalues')
        if strcmp(last_kind, 'alpha_virt') || strcmp(last_kind, 'beta_virt')
            blocks = local_push_block(blocks, cur_alpha_occ, cur_alpha_virt, cur_beta_occ, cur_beta_virt);
            cur_alpha_occ = []; cur_alpha_virt = []; cur_beta_occ = []; cur_beta_virt = [];
        end
        cur_alpha_occ = [cur_alpha_occ, local_extract_nums(ln)]; %#ok<AGROW>
        last_kind = 'alpha_occ';
    elseif contains(ln, 'Alpha virt. eigenvalues')
        cur_alpha_virt = [cur_alpha_virt, local_extract_nums(ln)]; %#ok<AGROW>
        last_kind = 'alpha_virt';
    elseif contains(ln, 'Beta  occ. eigenvalues') || contains(ln, 'Beta occ. eigenvalues')
        cur_beta_occ = [cur_beta_occ, local_extract_nums(ln)]; %#ok<AGROW>
        last_kind = 'beta_occ';
    elseif contains(ln, 'Beta virt. eigenvalues')
        cur_beta_virt = [cur_beta_virt, local_extract_nums(ln)]; %#ok<AGROW>
        last_kind = 'beta_virt';
    end
end
blocks = local_push_block(blocks, cur_alpha_occ, cur_alpha_virt, cur_beta_occ, cur_beta_virt);

if isempty(blocks)
    error('G16_orbital_energies: no orbital eigenvalue block found in %s', filename);
end

% -------------------------------------------------------------------------
% Select the requested step
% -------------------------------------------------------------------------
Nsteps = numel(blocks);
if ischar(step_req)
    if strcmpi(step_req, 'last'),  si = Nsteps;
    elseif strcmpi(step_req, 'first'), si = 1;
    else, error('G16_orbital_energies: step must be ''first'', ''last'', or an integer.');
    end
else
    si = round(step_req);
    if si < 1 || si > Nsteps
        error('G16_orbital_energies: step %d out of range [1, %d].', si, Nsteps);
    end
end

blk = blocks(si);
if isempty(blk.alpha_occ) || isempty(blk.alpha_virt)
    error('G16_orbital_energies: incomplete eigenvalue block (step %d) in %s', si, filename);
end

% -------------------------------------------------------------------------
% HOMO / LUMO
% -------------------------------------------------------------------------
HOMO_a = blk.alpha_occ(end);
LUMO_a = blk.alpha_virt(1);
has_beta = ~isempty(blk.beta_occ) && ~isempty(blk.beta_virt);

if has_beta
    HOMO_b = blk.beta_occ(end);
    LUMO_b = blk.beta_virt(1);
    HOMO = max(HOMO_a, HOMO_b);
    LUMO = min(LUMO_a, LUMO_b);
else
    HOMO_b = []; LUMO_b = [];
    HOMO = HOMO_a;
    LUMO = LUMO_a;
end
gap = LUMO - HOMO;

ha2eV = 27.211386245988;

% -------------------------------------------------------------------------
% Build output struct
% -------------------------------------------------------------------------
oe.alpha_occ  = blk.alpha_occ(:);
oe.alpha_virt = blk.alpha_virt(:);
oe.beta_occ   = blk.beta_occ(:);
oe.beta_virt  = blk.beta_virt(:);
oe.has_beta   = has_beta;

oe.HOMO       = HOMO;
oe.LUMO       = LUMO;
oe.gap        = gap;
oe.HOMO_alpha = HOMO_a;
oe.LUMO_alpha = LUMO_a;
oe.HOMO_beta  = HOMO_b;
oe.LUMO_beta  = LUMO_b;

oe.HOMO_eV = HOMO * ha2eV;
oe.LUMO_eV = LUMO * ha2eV;
oe.gap_eV  = gap  * ha2eV;

oe.step     = si;
oe.Nsteps   = Nsteps;
oe.filename = filename;

% -------------------------------------------------------------------------
% Print summary
% -------------------------------------------------------------------------
fprintf('\n── G16_orbital_energies (step %d/%d): %s ──\n', si, Nsteps, filename);
fprintf('  HOMO = %+.6f Ha  (%+.4f eV)\n', HOMO, oe.HOMO_eV);
fprintf('  LUMO = %+.6f Ha  (%+.4f eV)\n', LUMO, oe.LUMO_eV);
fprintf('  Gap  =  %.6f Ha  ( %.4f eV)\n', gap, oe.gap_eV);
if has_beta
    fprintf('  (open-shell: alpha HOMO/LUMO = %+.6f / %+.6f Ha, beta HOMO/LUMO = %+.6f / %+.6f Ha)\n', ...
        HOMO_a, LUMO_a, HOMO_b, LUMO_b);
end
fprintf('\n');

end % G16_orbital_energies


% =========================================================================
function blocks = local_push_block(blocks, a_occ, a_virt, b_occ, b_virt)
%LOCAL_PUSH_BLOCK  Appends an eigenvalue block if any data was collected.
if ~isempty(a_occ) || ~isempty(a_virt)
    blocks(end+1) = struct('alpha_occ', a_occ, 'alpha_virt', a_virt, ...
                            'beta_occ',  b_occ, 'beta_virt',  b_virt); %#ok<AGROW>
end
end

% -------------------------------------------------------------------------
function vals = local_extract_nums(ln)
%LOCAL_EXTRACT_NUMS  Pulls all signed decimal numbers out of a line, e.g.
%   " Alpha  occ. eigenvalues --  -19.12220 -19.12219 -14.36911"
%   -> [-19.12220 -19.12219 -14.36911]
vals = str2double(regexp(ln, '-?\d+\.\d+', 'match'));
end
