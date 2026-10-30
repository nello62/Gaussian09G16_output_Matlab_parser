function [angTable, ang] = G16_get_angle(mol, atoms, varargin)
% G16_GET_ANGLE  Bond angles i-j-k (degrees) from a molecular structure.
%
%   [angTable, ang] = G16_GET_ANGLE(mol, atoms)
%   [angTable, ang] = G16_GET_ANGLE(mol, atoms, 'Name', Value, ...)
%
%   INPUT
%     mol    - struct with fields .xyz [N x 3] (Angstrom) and, optionally,
%              .symbols {N x 1} (as returned by G16_structure,
%              G16_fchk_read(...).mol, or G_read_input), or a numeric
%              [N x 3] coordinate matrix
%     atoms  - [M x 3] atom indices (Gaussian numbering, 1-based); row m
%              defines the angle atoms(m,1)-atoms(m,2)-atoms(m,3), with
%              the vertex in the middle
%
%   OPTIONAL NAME-VALUE PAIRS
%     'SaveAs'   .xlsx or .csv file where the table is written
%     'Verbose'  true | false (default)  print the table
%
%   OUTPUT
%     angTable - table with columns Atom1, Sym1, Atom2, Sym2, Atom3, Sym3,
%                Angle_deg (symbols are '' if mol has no .symbols)
%     ang      - [M x 1] angles in degrees (0..180)
%
%   The angle is computed as atan2(|u x v|, u . v), with u = r1 - r2 and
%   v = r3 - r2, which is accurate also for angles close to 0 or 180 deg
%   (where acos loses precision).
%
%   EXAMPLES
%     mol = G16_structure('trans_optfreq.log');
%     G16_get_angle(mol, [4 7 8], 'Verbose', true)        % C4-N7=N8
%     T = G16_get_angle(mol, [4 7 8; 7 8 9]);              % both CNN angles
%
%   See also: G16_GET_DIHEDRAL, G16_GET_BOND_LENGTH, G16_STRUCTURE
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

p = inputParser;
addRequired(p,  'mol');
addRequired(p,  'atoms', @(x) isnumeric(x) && size(x, 2) == 3);
addParameter(p, 'SaveAs',  '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Verbose', false, @(x) islogical(x) || isnumeric(x));
parse(p, mol, atoms, varargin{:});
opt = p.Results;

[xyz, sym] = local_coords(mol, 'G16_get_angle');
check_indices(atoms, size(xyz, 1), 'G16_get_angle');

M   = size(atoms, 1);
ang = zeros(M, 1);
for m = 1:M
    u = xyz(atoms(m,1),:) - xyz(atoms(m,2),:);
    v = xyz(atoms(m,3),:) - xyz(atoms(m,2),:);
    ang(m) = atan2d(norm(cross(u, v)), dot(u, v));
end

a = atoms;
angTable = table(a(:,1), sym(a(:,1)), a(:,2), sym(a(:,2)), a(:,3), sym(a(:,3)), ang, ...
    'VariableNames', {'Atom1','Sym1','Atom2','Sym2','Atom3','Sym3','Angle_deg'});

if opt.Verbose
    fprintf('\n%-16s %10s\n', 'angle', 'deg');
    for m = 1:M
        lab = sprintf('%s%d-%s%d-%s%d', sym{a(m,1)}, a(m,1), sym{a(m,2)}, a(m,2), ...
                      sym{a(m,3)}, a(m,3));
        fprintf('%-16s %10.3f\n', lab, ang(m));
    end
    fprintf('\n');
end
if ~isempty(opt.SaveAs)
    writetable(angTable, char(opt.SaveAs));
    fprintf('Table saved to: %s\n', char(opt.SaveAs));
end
end

% -------------------------------------------------------------------------
function [xyz, sym] = local_coords(mol, fname)
if isnumeric(mol)
    xyz = mol;
    sym = repmat({''}, size(xyz, 1), 1);
elseif isstruct(mol) && isfield(mol, 'xyz')
    xyz = mol.xyz;
    if isfield(mol, 'symbols') && numel(mol.symbols) == size(xyz, 1)
        sym = cellstr(mol.symbols(:));
    else
        sym = repmat({''}, size(xyz, 1), 1);
    end
else
    error('%s: mol must be a struct with field .xyz or an N x 3 matrix.', fname);
end
if size(xyz, 2) ~= 3
    error('%s: coordinates must be an N x 3 matrix.', fname);
end
end

function check_indices(atoms, N, fname)
if any(atoms(:) < 1) || any(atoms(:) > N) || any(atoms(:) ~= round(atoms(:)))
    error('%s: atom indices must be integers between 1 and %d.', fname, N);
end
end
