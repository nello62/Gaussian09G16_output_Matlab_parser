function bt = G16_nbo_bonds(filename, varargin)
% G16_NBO_BONDS  Determines bond order (single/double/triple) from a
%                Gaussian NBO (Natural Bond Orbital) analysis.
%
%   bt = G16_NBO_BONDS(filename)
%   bt = G16_NBO_BONDS(filename, 'section', 'last')   % default
%   bt = G16_NBO_BONDS(filename, 'section', 'first')
%
%   Requires the source file to have been computed with the 'pop=nbo'
%   Gaussian keyword (or similar, e.g. 'pop=(nbo,savenbo)'); without it,
%   the output file simply does not contain NBO data and this function
%   raises an error -- bond order cannot be recovered after the fact
%   from geometry/energy alone, the NBO analysis has to have actually
%   run as part of the job.
%
%   Unlike G16_GET_BOND_LENGTH (a purely geometric covalent-radius
%   criterion) and G16_DRAW_MOLECULE's bond-order heuristic (a bond
%   length threshold for C-C/C-O pairs only, explicitly not derived from
%   any real Gaussian bond-order analysis), this reads the actual NBO
%   "Natural Bond Orbitals (Summary)" table and counts, for every atom
%   pair, how many bonding (BD) natural bond orbitals NBO assigned to
%   it: one BD orbital = single bond (sigma only), two = double bond
%   (sigma + pi), three = triple bond (sigma + 2 pi). This reflects
%   Gaussian's own NBO analysis, not a geometric guess.
%
%   Note: this reads the "Natural Bond Orbitals (Summary)" table that
%   'pop=nbo' always prints, not the separate "Wiberg bond index matrix"
%   (which additionally requires 'pop=(nbo,bndidx)' and is not parsed by
%   this function).
%
%   OUTPUT  table bt with columns:
%       Atom1, Atom2   int       atom indices (1-based, Atom1 < Atom2)
%       Sym1, Sym2     string    atomic symbols matching Atom1/Atom2
%       BondOrder      int       count of BD orbitals for this pair
%       Occupancy      double    summed NBO occupancy (~2 per BD orbital)
%
%   Optional parameters (Name-Value):
%       'section'  - 'last' (default) | 'first' -- which "Natural Bond
%                    Orbitals (Summary)" block to read, for files with
%                    more than one (e.g. multi-step opt+freq+polar jobs
%                    that repeat the NBO analysis at each step)
%       'Lines'    - pre-read cell array of file lines, to skip
%                    re-reading the file when it has already been read
%                    elsewhere. Default {} (read the file normally).
%
%   Example:
%       bt = G16_nbo_bonds('molecule_nbo.out');
%
%   See also G16_GET_BOND_LENGTH, G16_DRAW_MOLECULE.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

p = inputParser;
addRequired(p,  'filename', @ischar);
addParameter(p, 'section',  'last', @(x) ischar(x) && any(strcmpi(x, {'first','last'})));
addParameter(p, 'Lines',    {},     @iscell);
parse(p, filename, varargin{:});
section_req = lower(p.Results.section);

lines = p.Results.Lines;
if isempty(lines)
    if ~isfile(filename)
        error('G16_nbo_bonds: file not found: %s', filename);
    end
    fid = fopen(filename, 'r');
    raw = fread(fid, '*char')';
    fclose(fid);
    lines = strsplit(raw, newline, 'CollapseDelimiters', false);
    G16_check_gaussian_match(lines, filename);
end
N = numel(lines);

hdr_idx = find(~cellfun(@isempty, regexp(lines, 'Natural Bond Orbitals \(Summary\)', 'once')));
if isempty(hdr_idx)
    error('G16_nbo_bonds:notFound', ...
        ['No NBO analysis found in %s. This requires the source Gaussian ' ...
         'job to have been run with the ''pop=nbo'' keyword (or similar) ' ...
         '-- without it, bond-order data is not present in the file.'], filename);
end

if strcmp(section_req, 'last')
    k0 = hdr_idx(end);
else
    k0 = hdr_idx(1);
end

% Bonding (BD) NBO lines look like:
%     1. BD (   1) C   1 - H   2          1.99909    -0.50504
% Antibonding (BD*) lines are deliberately excluded: '\(' is required to
% follow 'BD' with only whitespace in between, so 'BD*(' cannot match.
bd_re = '^\s*\d+\.\s+BD\s*\(\s*\d+\)\s+([A-Za-z]+)\s*(\d+)\s*-\s*([A-Za-z]+)\s*(\d+)\s+([\d.]+)';

pairs = zeros(0, 2);
syms  = cell(0, 2);
occ   = zeros(0, 1);

for k = k0:N
    if ~isempty(regexp(lines{k}, 'Total Lewis', 'once'))
        break
    end
    tok = regexp(lines{k}, bd_re, 'tokens', 'once');
    if isempty(tok)
        continue
    end
    raw_a1 = str2double(tok{2});
    raw_a2 = str2double(tok{4});
    raw_s1 = tok{1};
    raw_s2 = tok{3};

    a1 = min(raw_a1, raw_a2);
    a2 = max(raw_a1, raw_a2);
    if raw_a1 == a1
        s1 = raw_s1; s2 = raw_s2;
    else
        s1 = raw_s2; s2 = raw_s1;
    end

    pairs(end+1, :) = [a1, a2]; %#ok<AGROW>
    syms(end+1, :)  = {s1, s2}; %#ok<AGROW>
    occ(end+1, 1)   = str2double(tok{5}); %#ok<AGROW>
end

if isempty(pairs)
    error('G16_nbo_bonds:noBonds', ...
        'NBO summary section found in %s but no bonding (BD) orbitals were read.', filename);
end

[uniq_pairs, first_idx, ic] = unique(pairs, 'rows', 'stable');
Nu = size(uniq_pairs, 1);

BondOrder = accumarray(ic, 1);
Occupancy = accumarray(ic, occ);

Atom1 = uniq_pairs(:, 1);
Atom2 = uniq_pairs(:, 2);
Sym1  = strings(Nu, 1);
Sym2  = strings(Nu, 1);
for u = 1:Nu
    Sym1(u) = syms{first_idx(u), 1};
    Sym2(u) = syms{first_idx(u), 2};
end

bt = table(Atom1, Sym1, Atom2, Sym2, BondOrder, Occupancy);
bt = sortrows(bt, {'Atom1', 'Atom2'});

fprintf('\n── G16_nbo_bonds: %s ──\n', filename);
fprintf('  %d bonds found from NBO analysis (section: %s)\n', Nu, section_req);
for u = 1:height(bt)
    fprintf('  %s%d - %s%d : order %d\n', bt.Sym1(u), bt.Atom1(u), bt.Sym2(u), bt.Atom2(u), bt.BondOrder(u));
end
fprintf('\n');

end % G16_nbo_bonds
