%% G16 Toolbox — Example #1: structure, vibrational mode, and charges
%
%   Demonstrates three core G16_*.m functions on a real DFT calculation:
%     1) G16_structure  + G16_draw_molecule  — 3D structure rendering
%     2) G16_nmodes     + G16_draw_mode      — a single vibrational mode
%     3) G16_charges                         — Mulliken charges + dipole
%
%   Data file: test_2.out — DFT geometry optimisation of Violacein.
%   Reference: G. Cassone et al., Phys. Chem. Chem. Phys.,
%              doi: 10.1039/d6cp01164k
%
%   Requirements: G16/ on the MATLAB path, test_2.out in the current
%   folder. Running this script creates test_2_mode_91.pdf in the
%   current folder.

clear
close all
clc

filename = 'test_2.out';

% Camera orientation shared by all three figures below, so the molecule
% is viewed from the same angle in every plot.
view_az = -16.4417;
view_el =  44.0197;

%% 1) Structure — load the optimised geometry and render it in 3D
mol = G16_structure(filename);

% Create the axes up front so it can be tweaked (view angle, export)
% after G16_draw_molecule has finished drawing into it.
ax1 = axes;
set(ax1, 'Visible', 'off');
G16_draw_molecule(mol, 'Title', 'Violacein', 'Ax', ax1);
set(ax1, 'View', [view_az, view_el]);

% Uncomment to export (ContentType 'vector' gives a crisp publication
% figure but is slow for a molecule this size; use 'image' for a quick
% low-effort preview instead).
% exportgraphics(ax1, 'test_2_mol.pdf', 'ContentType', 'vector');

%% 2) Vibrational mode — mode #91 is the main C=C stretch
nm = G16_nmodes(filename);

ax2 = axes;
set(ax2, 'Visible', 'off');
G16_draw_mode(mol, nm, 91);
set(ax2, 'View', [view_az, view_el]);

% Replace only the molecule-name line of the auto-generated two-line
% title; the second line (mode number, frequency, IR/Raman intensity)
% is left untouched.
t = get(gca, 'Title');
t.String(1) = {'Violacein'};
set(gca, 'Title', t);

exportgraphics(gca, 'test_2_mode_91.pdf', 'ContentType', 'vector');

%% 3) Atomic charges — Mulliken charges with a dipole moment overlay
G16_charges(filename, 'Plot', true, 'ShowDipole', true);
set(gca, 'View', [view_az, view_el]);

t = get(gca, 'Title');
t.String = {'Violacein - Mulliken Charges (atom)'};
set(gca, 'Title', t);

fprintf('\nDone. Figures: structure, mode 91, Mulliken charges.\n');
fprintf('Saved: test_2_mode_91.pdf\n');
