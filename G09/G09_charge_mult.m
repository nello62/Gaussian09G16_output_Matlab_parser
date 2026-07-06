function [charge, mult] = G09_charge_mult(filename)
% G09_CHARGE_MULT  Extracts molecular charge and spin multiplicity
%                  from a Gaussian 09 output file.
%
%   [charge, mult] = G09_CHARGE_MULT(filename)
%
%   Format is identical to G16:
%     "Charge =  0 Multiplicity = 1"
%
%   Example:
%       [q, m] = G09_charge_mult('indaco.log')   % q=0, m=1

if ~isfile(filename)
    error('G09_charge_mult: file not found: %s', filename);
end

lines = G09_read_lines(filename);

charge = NaN;
mult   = NaN;

for k = 1:numel(lines)
    ln = lines{k};
    tok = regexp(ln, 'Charge\s*=\s*([-\d]+)\s+Multiplicity\s*=\s*(\d+)', ...
                 'tokens', 'once');
    if ~isempty(tok)
        charge = str2double(tok{1});
        mult   = str2double(tok{2});
        break
    end
end

if isnan(charge)
    error('G09_charge_mult: "Charge = ... Multiplicity = ..." not found in %s', filename);
end

fprintf('Charge = %+d   Multiplicity = %d\n', charge, mult);

end  % G09_charge_mult
