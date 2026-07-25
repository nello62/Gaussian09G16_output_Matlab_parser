%% G16 Toolbox — Example #3: Energetics
%
%   Demonstrates G16_orbital_energies, G16_energy, and G16_convergence
%   on a real DFT calculation, then combines their results into a
%   one-line summary — a small taste of what G16_read_all/
%   G16_write_report do at full scale across every extraction function.
%   Also renders the HOMO-LUMO diagram via G16_draw_orbital.
%
%   Data file: test_2.out — DFT geometry optimisation of Violacein.
%   Reference: G. Cassone et al., Phys. Chem. Chem. Phys.,
%              doi: 10.1039/d6cp01164k
%
%   Requirements: G16/ on the MATLAB path, test_2.out in the current
%   folder. Each function already prints its own formatted summary to
%   the command window (shown below as a reference for what to expect
%   on this file); this script does not repeat those printouts. Running
%   it creates test_2_orbital_diagram.pdf in the current folder.

clear; close all; clc
filename = 'test_2.out';
%% 1) HOMO-LUMO gap — detailed data in the oe structure
oe = G16_orbital_energies(filename);
% Console output for this file:
%   ── G16_orbital_energies (step 3/3): test_2.out ──
%     HOMO = -0.190210 Ha  (-5.1759 eV)
%     LUMO = -0.096890 Ha  (-2.6365 eV)
%     Gap  =  0.093320 Ha  ( 2.5394 eV)

% Orbital energy-level diagram (title defaults to the source filename)
G16_draw_orbital(oe);
exportgraphics(gcf, 'test_2_orbital_diagram.pdf', 'ContentType', 'vector');
%% 2) Energy and thermochemistry — detailed data in the en structure
en = G16_energy(filename);
% Console output for this file:
%   ── G16_energy: test_2.out ──
%     Method : RB3LYP
%     SCF    : -1160.20276011  Ha
%     ZPE    : +0.29182000  Ha   (766.17 kJ/mol)
%     E0+ZPE : -1159.91094000  Ha
%     H      : -1159.88988100  Ha
%     G      : -1159.96045200  Ha
%     S      :  621.4460  J/(mol·K)
%     T = 298.15 K,  P = 1.00000 atm

%% 3) Geometry-optimisation convergence — detailed data in the cv structure
cv = G16_convergence(filename);
% Console output for this file:
%   G16_convergence: test_2.out
%     Steps read  : 11
%     Converged   : YES  (step 10)
%     Last step   : MaxF=1.00e-06 (thr 1.50e-05)  RMSF=0.00e+00 (thr 1.00e-05)
%                   MaxD=1.24e-04 (thr 6.00e-05)  RMSD=2.50e-05 (thr 4.00e-05)

%% 4) Combine all three into a one-line summary
if cv.converged
    conv_label = sprintf('YES (step %d/%d)', cv.conv_step, cv.Nsteps);
else
    conv_label = sprintf('NO (%d steps read)', cv.Nsteps);
end
fprintf('\nSummary for %s:\n', filename);
fprintf('  HOMO-LUMO gap : %.3f eV\n', oe.gap_eV);
if en.has_thermo
    fprintf('  Gibbs energy  : %.2f kJ/mol\n', en.G_kJ);
else
    fprintf('  Gibbs energy  : n/a (no frequency/thermochemistry job)\n');
end
fprintf('  Converged     : %s\n', conv_label);
