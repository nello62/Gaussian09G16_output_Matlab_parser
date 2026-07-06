function T = G16_list()
% G16_LIST  Lists all G16_*.m functions in this toolbox with a one-line description.
%
%   G16_LIST() prints a formatted list of every G16_*.m file found in the
%   same folder as this function (name + H1 description line), sorted
%   alphabetically.
%
%   T = G16_LIST() also returns the list as a table with columns:
%       .Name          function name (without .m)
%       .Description   H1 comment line (first comment line after the
%                       function declaration), '' if not found
%       .File          full path to the .m file
%
%   Example:
%       G16_list();
%       T = G16_list();
%       T(contains(T.Description, 'dipole', 'IgnoreCase', true), :)

thisFile = mfilename('fullpath');
folder   = fileparts(thisFile);

files = dir(fullfile(folder, 'G16_*.m'));
names = sort({files.name});

Name        = strings(numel(names), 1);
Description = strings(numel(names), 1);
File        = strings(numel(names), 1);

for k = 1:numel(names)
    fpath = fullfile(folder, names{k});
    [~, fnameNoExt] = fileparts(names{k});
    Name(k)        = string(fnameNoExt);
    Description(k) = string(local_h1_line(fpath));
    File(k)        = string(fpath);
end

T = table(Name, Description, File);

fprintf('\n── G16 Toolbox — %d function(s) in %s ──\n\n', numel(names), folder);
nameW = max(strlength(Name));
for k = 1:numel(names)
    fprintf('  %-*s  %s\n', nameW, Name(k), Description(k));
end
fprintf('\n');

end % G16_list


% =========================================================================
function h1 = local_h1_line(fpath)
%LOCAL_H1_LINE  Returns the H1 comment line (first % line right after the
%   function declaration), stripped of the leading '%'. '' if not found.
h1 = '';
fid = fopen(fpath, 'r');
if fid == -1
    return
end
firstLine = true;
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    ln = strtrim(ln);
    if firstLine
        firstLine = false;
        continue   % skip the "function ... = ..." declaration line
    end
    if isempty(ln) || ln(1) ~= '%'
        break
    end
    h1 = strtrim(ln(2:end));
    break
end
fclose(fid);
end
