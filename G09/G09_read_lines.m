function lines = G09_read_lines(filename)
% G09_READ_LINES  Read a Gaussian 09 output file and return a cell array of lines.
%
%   lines = G09_READ_LINES(filename)
%
%   Handles:
%     - CRLF line endings (Windows G09W output)
%     - latin-1 (ISO-8859-1) encoding
%
%   Returns a cell array of strings with line endings stripped.

if ~isfile(filename)
    error('G09_read_lines: file not found: %s', filename);
end

fid  = fopen(filename, 'r', 'n', 'ISO-8859-1');
raw  = fread(fid, '*char')';
fclose(fid);

% Normalise line endings
raw   = strrep(raw, sprintf('\r\n'), newline);
raw   = strrep(raw, sprintf('\r'),   newline);
lines = strsplit(raw, newline);

G09_check_gaussian_match(lines, filename);
end
