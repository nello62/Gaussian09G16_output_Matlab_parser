function gv = G09_gaussian_version(filename)
% G09_GAUSSIAN_VERSION  Detects which Gaussian version/revision produced a
%                       .out/.log/.fchk file.
%
%   gv = G09_GAUSSIAN_VERSION(filename)
%
%   For .out/.log files this reads the "Gaussian NN, Revision X.YY,"
%   citation line that Gaussian prints near the top of every output file
%   (works for any NN: 09, 16, ...).
%
%   .fchk files do NOT store this information at all — the formatted
%   checkpoint format has no version/provenance field. In that case this
%   function looks for a sibling .log/.out file with the same base name
%   in the same folder (as typically produced alongside the .fchk by the
%   same job) and reads the version from there instead. If no sibling
%   file is found, .major/.revision/.full come back empty and a warning
%   is issued.
%
%   OUTPUT  struct gv with fields:
%       .major      double   Gaussian major version, e.g. 9, 16 ([] if unknown)
%       .revision   char     revision string, e.g. 'A.02', 'C.01' ('' if unknown)
%       .full       char     citation line as printed, e.g.
%                            'Gaussian 09, Revision A.02' ('' if unknown)
%       .source     char     'out/log' | 'fchk-sibling:<path>' | 'unknown'
%       .filename   char     the input filename
%
%   Example:
%       gv = G09_gaussian_version('a1.out');       % gv.major = 16, gv.revision = 'C.01'
%       gv = G09_gaussian_version('molecule.fchk'); % falls back to molecule.out/.log if present

if ~isfile(filename)
    error('G09_gaussian_version: file not found: %s', filename);
end

[~, ~, ext] = fileparts(filename);

if strcmpi(ext, '.fchk')
    [major, revision, full] = local_read_version(filename);
    if ~isempty(major)
        source = 'out/log';   % unexpected for a genuine .fchk, but kept for completeness
    else
        sibling = local_find_sibling(filename);
        if ~isempty(sibling)
            [major, revision, full] = local_read_version(sibling);
            if ~isempty(major)
                source = sprintf('fchk-sibling:%s', sibling);
            else
                source = 'unknown';
            end
        else
            source = 'unknown';
        end
    end
else
    [major, revision, full] = local_read_version(filename);
    if ~isempty(major)
        source = 'out/log';
    else
        source = 'unknown';
    end
end

if isempty(major)
    warning('G09_gaussian_version:notFound', ...
        ['Could not determine the Gaussian version for %s ', ...
         '(no "Gaussian NN, Revision ..." line found, and no readable ', ...
         'sibling .log/.out file).'], filename);
end

gv.major    = major;
gv.revision = revision;
gv.full     = full;
gv.source   = source;
gv.filename = filename;

fprintf('%s -> %s\n', filename, local_display_string(gv));

end % G09_gaussian_version


% =========================================================================
function [major, revision, full] = local_read_version(fpath)
%LOCAL_READ_VERSION  Scans the first lines of FPATH for the "Gaussian NN,
%   Revision X.YY," citation line. Returns [] / '' for all outputs if not
%   found within the scan window (the line always appears very early in a
%   genuine Gaussian output file).
major = []; revision = ''; full = '';
fid = fopen(fpath, 'r');
if fid == -1
    return
end
n = 0;
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    n = n + 1;
    tok = regexp(ln, 'Gaussian\s+(\d+)\s*,\s*Revision\s+([\w.]+)', 'tokens', 'once');
    if ~isempty(tok)
        major    = str2double(tok{1});
        revision = tok{2};
        full     = regexprep(strtrim(ln), ',\s*$', '');
        break
    end
    if n > 500   % the citation line is always within the first few dozen lines
        break
    end
end
fclose(fid);
end

% -------------------------------------------------------------------------
function sibling = local_find_sibling(fchk_path)
%LOCAL_FIND_SIBLING  Looks for a .log/.out file with the same base name as
%   FCHK_PATH in the same folder. Returns '' if none exists.
sibling = '';
[folder, base] = fileparts(fchk_path);
candidates = {'.log', '.out', '.LOG', '.OUT'};
for i = 1:numel(candidates)
    cand = fullfile(folder, [base, candidates{i}]);
    if isfile(cand)
        sibling = cand;
        return
    end
end
end

% -------------------------------------------------------------------------
function s = local_display_string(gv)
if isempty(gv.major)
    s = 'unknown Gaussian version';
else
    s = sprintf('%s (revision %s)', gv.full, gv.revision);
    if ~strcmp(gv.source, 'out/log')
        s = sprintf('%s [from %s]', s, gv.source);
    end
end
end
