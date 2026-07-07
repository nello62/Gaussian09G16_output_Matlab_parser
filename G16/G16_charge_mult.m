function [charge, mult] = G16_charge_mult(filename, varargin)
% G16_CHARGE_MULT  Extracts charge and spin multiplicity from a .out/.log file
%                  .out/.log file.
%
%   [charge, mult] = G16_CHARGE_MULT(filename)
%
%   OUTPUT:
%       charge  int   molecular charge (e.g. 0, -1, +2, ...)
%       mult    int   spin multiplicity (1=singlet, 2=doublet, 3=triplet, ...)
%
%   Optional parameters:
%       'Lines'  - pre-read cell array of file lines, to skip re-reading
%                  the file when it has already been read elsewhere (e.g.
%                  G16_READ_ALL). Default {} (read the file normally).
%
%   Example:
%       [q, m] = G16_charge_mult('V_E00t.out')
%       % q = 0,  m = 1

p = inputParser;
addRequired(p,  'filename', @ischar);
addParameter(p, 'Lines',    {}, @iscell);
parse(p, filename, varargin{:});

lines = p.Results.Lines;
if isempty(lines)
    if ~isfile(filename)
        error('G16_charge_mult: file not found: %s', filename);
    end
    fid  = fopen(filename, 'r');
    raw  = fread(fid, '*char')';
    fclose(fid);
    lines = strsplit(raw, newline);
end

charge = NaN;
mult   = NaN;

for k = 1:numel(lines)
    ln = lines{k};
    % "Charge =  0 Multiplicity = 1"
    tok = regexp(ln, 'Charge\s*=\s*([-\d]+)\s+Multiplicity\s*=\s*(\d+)', ...
                 'tokens', 'once');
    if ~isempty(tok)
        charge = str2double(tok{1});
        mult   = str2double(tok{2});
        break   % use first occurrence (from the input section)
    end
end

if isnan(charge)
    error('G16_charge_mult: "Charge = ... Multiplicity = ..." line not found in %s', filename);
end

fprintf('Charge = %+d   Multiplicity = %d\n', charge, mult);

end
