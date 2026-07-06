function bondTable = G16_get_bond_length(mol, varargin)
%G16_GET_BOND_LENGTH  Costruisce la tabella delle distanze di legame di una molecola.
%
%   bondTable = G16_GET_BOND_LENGTH(mol) analizza la struttura MOL (come
%   quelle prodotte dal toolbox G16/G16, con campi .symbols [1xN cell o
%   string array] e .xyz [Nx3 double]) e restituisce una TABLE con tutti
%   i legami individuati tramite un criterio basato sui raggi covalenti.
%
%   bondTable = G16_GET_BOND_LENGTH(mol, 'Name', Value, ...) consente di
%   specificare le seguenti opzioni:
%
%     'Tolerance'   fattore moltiplicativo sulla somma dei raggi
%                   covalenti per decidere se due atomi sono legati
%                   (default: 1.15)
%     'IncludeH'    true/false, se includere i legami che coinvolgono
%                   idrogeni (default: true)
%     'SortBy'      'distance' (default) oppure 'atom' per ordinare la
%                   tabella per distanza crescente o per indice atomico
%     'SaveAs'      percorso file (.xlsx o .csv) su cui salvare la
%                   tabella; se omesso non viene salvato nulla
%
%   La tabella in output ha le colonne:
%     Atom1, Sym1, Atom2, Sym2, Distance_Ang
%
%   Esempio:
%     load('M.mat');
%     T = G16_get_bond_length(molzn, 'SortBy', 'atom');
%     T = G16_get_bond_length(mol3tpy, 'IncludeH', false, ...
%                              'SaveAs', 'bonds_3tpy.xlsx');
%
%   Nota: il criterio di legame è puramente geometrico (raggi covalenti),
%   non deriva dalla connettività/topologia calcolata da Gaussian: va
%   quindi bene per un'analisi strutturale rapida, non per la matrice di
%   connettività "ufficiale" dell'output G16/G16.

p = inputParser;
addParameter(p, 'Tolerance', 1.15, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'IncludeH', true, @islogical);
addParameter(p, 'SortBy', 'distance', @(x) any(strcmpi(x, {'distance','atom'})));
addParameter(p, 'SaveAs', '', @ischar);
parse(p, varargin{:});
opt = p.Results;

symbols = cellstr(mol.symbols(:));
xyz     = mol.xyz;
Natoms  = numel(symbols);

if size(xyz,1) ~= Natoms
    error('G16_get_bond_length:sizeMismatch', ...
        'Il numero di simboli (%d) non corrisponde al numero di righe di xyz (%d).', ...
        Natoms, size(xyz,1));
end

R = covalent_radii();

Atom1 = zeros(0,1); Sym1 = strings(0,1);
Atom2 = zeros(0,1); Sym2 = strings(0,1);
Distance_Ang = zeros(0,1);

for i = 1:Natoms-1
    si = symbols{i};
    if ~opt.IncludeH && strcmpi(si,'H'), continue; end
    ri = lookup_radius(R, si);

    for j = i+1:Natoms
        sj = symbols{j};
        if ~opt.IncludeH && strcmpi(sj,'H'), continue; end
        rj = lookup_radius(R, sj);

        d = norm(xyz(i,:) - xyz(j,:));
        if d <= opt.Tolerance * (ri + rj)
            Atom1(end+1,1) = i;              %#ok<AGROW>
            Sym1(end+1,1)  = string(si);      %#ok<AGROW>
            Atom2(end+1,1) = j;               %#ok<AGROW>
            Sym2(end+1,1)  = string(sj);      %#ok<AGROW>
            Distance_Ang(end+1,1) = d;        %#ok<AGROW>
        end
    end
end

bondTable = table(Atom1, Sym1, Atom2, Sym2, Distance_Ang);

switch lower(opt.SortBy)
    case 'distance'
        bondTable = sortrows(bondTable, 'Distance_Ang');
    case 'atom'
        bondTable = sortrows(bondTable, {'Atom1','Atom2'});
end

if ~isempty(opt.SaveAs)
    writetable(bondTable, opt.SaveAs);
    fprintf('Tabella dei legami salvata in: %s\n', opt.SaveAs);
end

end % function G16_get_bond_length


% ------------------------------------------------------------------
function R = covalent_radii()
%COVALENT_RADII  Raggi covalenti (Angstrom) - Cordero et al., Dalton Trans. 2008
elements = {'H','C','N','O','F','S','P','Cl','Br','I', ...
            'Si','B','Na','K','Mg','Ca', ...
            'Au','Ag','Cu','Zn','Cd','Pt','Pd','Ni','Fe','Co'};
radii    = [0.31 0.76 0.71 0.66 0.57 1.05 1.07 1.02 1.20 1.39, ...
            1.11 0.84 1.66 2.03 1.41 1.76, ...
            1.36 1.45 1.32 1.22 1.44 1.36 1.39 1.24 1.32 1.26];
R = containers.Map(elements, num2cell(radii));
end

function r = lookup_radius(R, sym)
if isKey(R, sym)
    r = R(sym);
else
    warning('G16_get_bond_length:unknownElement', ...
        'Elemento "%s" non presente nella tabella dei raggi covalenti: uso raggio di default 1.50 A.', sym);
    r = 1.50;
end
end
