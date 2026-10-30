function [dihTable, phi] = G09_get_dihedral(mol, atoms, varargin)
% G09_GET_DIHEDRAL  Dihedral (torsion) angles i-j-k-l (degrees) from a molecular structure.
%
%   [dihTable, phi] = G09_GET_DIHEDRAL(mol, atoms)
%   [dihTable, phi] = G09_GET_DIHEDRAL(mol, atoms, 'Name', Value, ...)
%
%   INPUT
%     mol    - struct with fields .xyz [N x 3] (Angstrom) and, optionally,
%              .symbols {N x 1} (as returned by G09_structure,
%              G09_fchk_read(...).mol, or G_read_input), or a numeric
%              [N x 3] coordinate matrix
%     atoms  - [M x 4] atom indices (Gaussian numbering, 1-based); row m
%              defines the dihedral atoms(m,1)-(m,2)-(m,3)-(m,4), i.e. the
%              angle between the planes (1,2,3) and (2,3,4) seen along the
%              central bond 2-3
%
%   OPTIONAL NAME-VALUE PAIRS
%     'Signed'   true (default): signed dihedral in (-180, 180], IUPAC
%                convention (positive = clockwise rotation of the front
%                bond 1-2 onto the back bond 3-4 seen along 2->3);
%                false: absolute value in [0, 180]
%     'SaveAs'   .xlsx or .csv file where the table is written
%     'Verbose'  true | false (default)  print the table
%
%   OUTPUT
%     dihTable - table with columns Atom1, Sym1, ..., Atom4, Sym4,
%                Dihedral_deg (symbols are '' if mol has no .symbols)
%     phi      - [M x 1] dihedral angles in degrees
%
%   With b1 = r2 - r1, b2 = r3 - r2, b3 = r4 - r3, and the plane normals
%   n1 = b1 x b2, n2 = b2 x b3:
%       phi = atan2( (n1 x n2) . b2/|b2| , n1 . n2 )
%   The dihedral is undefined (NaN, with a warning) if three consecutive
%   atoms are collinear (e.g. at a linear CNN angle).
%
%   EXAMPLES
%     mol = G09_structure('cis_optfreq.log');
%     G09_get_dihedral(mol, [4 7 8 9], 'Verbose', true)    % CNNC
%     [~, phi] = G09_get_dihedral(mol, [4 7 8 9], 'Signed', false);
%
%   See also: G09_GET_ANGLE, G09_GET_BOND_LENGTH, G09_STRUCTURE
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

p = inputParser;
addRequired(p,  'mol');
addRequired(p,  'atoms', @(x) isnumeric(x) && size(x, 2) == 4);
addParameter(p, 'Signed',  true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SaveAs',  '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Verbose', false, @(x) islogical(x) || isnumeric(x));
parse(p, mol, atoms, varargin{:});
opt = p.Results;

[xyz, sym] = local_coords(mol, 'G09_get_dihedral');
check_indices(atoms, size(xyz, 1), 'G09_get_dihedral');

M   = size(atoms, 1);
phi = zeros(M, 1);
tol = 1e-8;
for m = 1:M
    r  = xyz(atoms(m,:), :);
    b1 = r(2,:) - r(1,:);  b2 = r(3,:) - r(2,:);  b3 = r(4,:) - r(3,:);
    n1 = cross(b1, b2);    n2 = cross(b2, b3);
    if norm(n1) < tol * norm(b1) * norm(b2) || norm(n2) < tol * norm(b2) * norm(b3)
        warning('G09_get_dihedral:collinear', ...
            'Dihedral %d-%d-%d-%d undefined: three consecutive atoms are collinear.', ...
            atoms(m,1), atoms(m,2), atoms(m,3), atoms(m,4));
        phi(m) = NaN;
        continue
    end
    phi(m) = atan2d(dot(cross(n1, n2), b2 / norm(b2)), dot(n1, n2));
end
if ~opt.Signed, phi = abs(phi); end

a = atoms;
dihTable = table(a(:,1), sym(a(:,1)), a(:,2), sym(a(:,2)), ...
                 a(:,3), sym(a(:,3)), a(:,4), sym(a(:,4)), phi, ...
    'VariableNames', {'Atom1','Sym1','Atom2','Sym2','Atom3','Sym3', ...
                      'Atom4','Sym4','Dihedral_deg'});

if opt.Verbose
    fprintf('\n%-22s %10s\n', 'dihedral', 'deg');
    for m = 1:M
        lab = sprintf('%s%d-%s%d-%s%d-%s%d', sym{a(m,1)}, a(m,1), sym{a(m,2)}, a(m,2), ...
                      sym{a(m,3)}, a(m,3), sym{a(m,4)}, a(m,4));
        fprintf('%-22s %10.3f\n', lab, phi(m));
    end
    fprintf('\n');
end
if ~isempty(opt.SaveAs)
    writetable(dihTable, char(opt.SaveAs));
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
