function aobasis = g_read_aobasis_from_fchk(filename)
%G_READ_AOBASIS_FROM_FCHK  Minimal, standalone .fchk section parser that
%   reads only the raw atomic-orbital basis-set sections (shell types,
%   primitive exponents, contraction coefficients, shell centres, ...)
%   needed to evaluate a molecular orbital or the electron density on a
%   real-space grid. Does not depend on, and does not duplicate the full
%   scope of, G09/G16_fchk_read -- kept deliberately self-contained so
%   this G_Utility function does not require any change to the core
%   toolbox. Shared by G_draw_mo_surface and G_draw_density_surface.
%
%   aobasis = g_read_aobasis_from_fchk(filename)
%
%   Output: struct with fields .shell_types, .n_prim_per_shell,
%   .shell_to_atom, .prim_exponents, .contraction_coeff,
%   .sp_contraction_coeff ([] if no SP shells), .shell_coords_bohr
%   ([Nshell x 3], [] if the section was absent/malformed).
%
%   .fchk format:
%     Scalar:  "<label padded>  I|R|C  <value>"
%     Array:   "<label padded>  I|R|C  N=  <count>" followed by data lines
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

fid   = fopen(filename, 'r');
raw   = fread(fid, '*char')';
fclose(fid);
lines = strsplit(raw, newline);
N     = numel(lines);

sec_names  = {};
sec_types  = {};
sec_counts = [];
sec_idx    = [];

header_re = '^(.{1,43}?)\s{1,3}(I|R|C)\s+(N=\s*(\d+)|[-\d.E+]+)';

for k = 1:N
    ln = lines{k};
    if numel(ln) < 45, continue; end
    tok = regexp(ln, header_re, 'tokens', 'once');
    if isempty(tok), continue; end
    name  = strtrim(tok{1});
    dtype = tok{2};
    rest  = strtrim(tok{3});
    if startsWith(rest, 'N=')
        nvals = str2double(strtrim(rest(3:end)));
    else
        nvals = 0;
    end
    sec_names{end+1}  = name;  %#ok<AGROW>
    sec_types{end+1}  = dtype; %#ok<AGROW>
    sec_counts(end+1) = nvals; %#ok<AGROW>
    sec_idx(end+1)    = k;     %#ok<AGROW>
end

    function val = read_sec(keyword)
        idx_match = find(~cellfun(@isempty, ...
            regexpi(sec_names, ['^', regexptranslate('wildcard', keyword), '$'])));
        if isempty(idx_match)
            idx_match = find(~cellfun(@isempty, ...
                cellfun(@(s) strfind(lower(s), lower(keyword)), sec_names, 'UniformOutput', false)));
        end
        if isempty(idx_match)
            val = [];
            return
        end
        mi    = idx_match(end);
        k0    = sec_idx(mi);
        dtype = sec_types{mi};
        nvals = sec_counts(mi);

        if nvals == 0
            tok = regexp(lines{k0}, '(I|R|C)\s+([-\d.E+]+)\s*$', 'tokens', 'once');
            if isempty(tok)
                val = [];
            else
                val = str2double(tok{2});
            end
            return
        end

        vals = [];
        k2   = k0 + 1;
        while numel(vals) < nvals && k2 <= N
            ln2 = strtrim(lines{k2});
            if isempty(ln2), k2 = k2+1; continue; end
            if numel(ln2) > 44 && ~isempty(regexp(ln2, header_re, 'once'))
                break
            end
            v = sscanf(ln2, '%f');
            vals = [vals; v]; %#ok<AGROW>
            k2 = k2 + 1;
        end
        val = vals(1:min(end, nvals));
    end

shell_types      = round(read_sec('Shell types'));
n_prim_per_shell = round(read_sec('Number of primitives per shell'));
shell_to_atom    = round(read_sec('Shell to atom map'));
prim_exponents    = read_sec('Primitive exponents');
contraction_coeff = read_sec('Contraction coefficients');
sp_contraction_coeff = read_sec('P(S=P) Contraction coefficients');   % [] if no SP shells
shell_coords_raw = read_sec('Coordinates of each shell');

aobasis.shell_types          = shell_types(:);
aobasis.n_prim_per_shell     = n_prim_per_shell(:);
aobasis.shell_to_atom        = shell_to_atom(:);
aobasis.prim_exponents       = prim_exponents(:);
aobasis.contraction_coeff    = contraction_coeff(:);
aobasis.sp_contraction_coeff = sp_contraction_coeff(:);
if numel(shell_coords_raw) == 3*numel(shell_types)
    aobasis.shell_coords_bohr = reshape(shell_coords_raw, 3, numel(shell_types))';   % [Nshell x 3]
else
    aobasis.shell_coords_bohr = [];
end

end % g_read_aobasis_from_fchk
