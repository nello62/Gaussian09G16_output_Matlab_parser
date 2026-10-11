function [cpkColors, covRadii] = G_element_table()
%G_ELEMENT_TABLE  Shared CPK colour and covalent-radius lookup tables.
%
%   [CPKCOLORS, COVRADII] = G_ELEMENT_TABLE() returns two containers.Map,
%   both keyed by proper-case element symbol (e.g. 'C', 'Cl'):
%       CPKCOLORS  - [R G B] triplet, CPK display colour convention
%       COVRADII   - covalent radius in Angstrom (Cordero et al., Dalton
%                    Trans. 2008)
%
%   The Gaussian input-file format does not differ between Gaussian 09 and
%   Gaussian 16, so (like G_READ_INPUT/G_SPECTRA_NM) this is a single
%   shared function rather than separate G09_/G16_ versions -- an
%   identical copy is kept in both the G09/ and G16/ folders.
%
%   Consolidates the two previously-independent, differently-sized
%   element tables embedded in G16_GET_BOND_LENGTH (26 elements, radii
%   only) and G16_DRAW_MOLECULE (21 elements, colour + radii) into one
%   shared source, so new code (currently: G_MOLECULE_BUILDER) has a
%   single table to build on instead of adding a third copy. The two
%   existing embedded tables are left exactly as they are for now (not
%   refactored to call this) -- see G_MOLECULE_BUILDER's own header
%   comment for why.
%
%   Callers are responsible for their own ISKEY(...) check and fallback
%   default on an unknown symbol (matching the existing
%   GET_RADIUS_LOCAL/LOOKUP_RADIUS precedent in G16_draw_molecule.m/
%   G16_get_bond_length.m) -- this function only returns the two Maps,
%   since neither existing fallback (1.50 vs 0.80 Angstrom default
%   radius) is more "right" out of context.
%
%   Example:
%       [cpk, rad] = G_element_table();
%       r = rad('C');           % 0.76
%       if isKey(cpk, 'Fe'), c = cpk('Fe'); else, c = [0.65 0.20 0.80]; end
%
%   See also G16_GET_BOND_LENGTH, G16_DRAW_MOLECULE, G_MOLECULE_BUILDER.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

    radiusElements = {'H','C','N','O','F','S','P','Cl','Br','I', ...
                      'Si','B','Na','K','Mg','Ca', ...
                      'Au','Ag','Cu','Zn','Cd','Pt','Pd','Ni','Fe','Co'};
    radii = [0.31 0.76 0.71 0.66 0.57 1.05 1.07 1.02 1.20 1.39, ...
             1.11 0.84 1.66 2.03 1.41 1.76, ...
             1.36 1.45 1.32 1.22 1.44 1.36 1.39 1.24 1.32 1.26];
    covRadii = containers.Map(radiusElements, num2cell(radii));

    colorElements = {'H',  'C',  'N',  'O',  'F',  'P',  'S',  'Cl', 'Br', 'I', ...
                     'Au', 'Ag', 'Fe', 'Zn', 'Ca', 'Mg', 'Na', 'K',  'Si', 'B', 'Cu'};
    colors = { ...
        [0.60 0.80 1.00], ...   % H  - light blue
        [0.30 0.30 0.30], ...   % C  - dark grey
        [0.10 0.30 0.90], ...   % N  - blue
        [0.90 0.10 0.10], ...   % O  - red
        [0.20 0.80 0.20], ...   % F  - green
        [1.00 0.50 0.00], ...   % P  - orange
        [1.00 0.85 0.00], ...   % S  - yellow
        [0.20 0.85 0.20], ...   % Cl - light green
        [0.55 0.20 0.10], ...   % Br - reddish brown
        [0.45 0.00 0.65], ...   % I  - purple
        [1.00 0.82 0.14], ...   % Au - gold
        [0.75 0.75 0.75], ...   % Ag - silver
        [0.80 0.40 0.00], ...   % Fe - orange-brown
        [0.50 0.70 0.50], ...   % Zn - teal green
        [0.60 0.60 0.60], ...   % Ca - grey
        [0.50 0.80 0.20], ...   % Mg - lime green
        [0.65 0.40 0.90], ...   % Na - light purple
        [0.55 0.20 0.85], ...   % K  - purple
        [0.60 0.50 0.40], ...   % Si - beige
        [1.00 0.65 0.50], ...   % B  - salmon
        [0.72 0.45 0.20]};      % Cu - copper
    cpkColors = containers.Map(colorElements, colors);
end
