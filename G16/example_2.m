%% G16 Toolbox — Example #2: Infrared and Raman spectra
%
%   Demonstrates G16_spectra on a real DFT calculation, two ways:
%     1) Manual plot — build a single Raman continuum plot directly
%        from the raw sp.x/sp.Raman_cont fields returned by G16_spectra.
%     2) Built-in plot — let G16_spectra draw its own two-panel
%        Raman+IR figure ('plot', true), then restyle the axes/labels
%        afterwards (larger, bold text) for a print-ready figure.
%
%   Data file: test_2.out — DFT geometry optimisation of Violacein.
%   Reference: G. Cassone et al., Phys. Chem. Chem. Phys.,
%              doi: 10.1039/d6cp01164k
%
%   Requirements: G16/ on the MATLAB path, test_2.out in the current
%   folder. Running this script creates two PDFs in the current
%   folder: test_2_raman_manual.pdf, test_2_ir_raman.pdf.

clear; close all; clc
filename = 'test_2.out';
%% 1) Manual plot — Raman continuum only, from the raw sp fields
sp = G16_spectra(filename);
fig1 = figure('Color', 'white');
plot(sp.x, sp.Raman_cont, 'LineWidth', 2);
xlabel('Wavenumber (cm^{-1})', 'FontWeight', 'bold');
ylabel('Raman activity (Å^4 AMU^{-1})', 'FontWeight', 'bold');
set(gca, 'FontSize', 12);
exportgraphics(fig1, 'test_2_raman_manual.pdf', 'ContentType', 'vector');
%% 2) Built-in two-panel plot (Raman + IR), then restyled
% 'FWHM' broadens the Lorentzian lines and 'xmin' crops the noisy
% low-frequency tail; 'plot', true triggers G16_spectra's own
% two-subplot Raman/IR figure (see G16_plot_spectra, local to
% G16_spectra.m).
G16_spectra(filename, 'FWHM', 15, 'xmin', 400, 'plot', true);
fig2 = gcf;
% G16_spectra creates the Raman axes first, then the IR axes below it;
% findobj returns them in reverse creation order, so flip to restore
% [Raman; IR].
axs      = flipud(findobj(fig2, 'Type', 'axes'));
ax_raman = axs(1);
ax_ir    = axs(2);
% Upsize the auto-generated labels (FontSize 10, regular weight) to a
% bolder, larger style for a print-ready figure.
set([ax_raman, ax_ir], 'FontSize', 12);
xlabel(ax_raman, 'Wavenumber (cm^{-1})',         'FontWeight', 'bold', 'FontSize', 12);
ylabel(ax_raman, 'Raman activity (Å^4 AMU^{-1})', 'FontWeight', 'bold', 'FontSize', 12);
xlabel(ax_ir,    'Wavenumber (cm^{-1})',         'FontWeight', 'bold', 'FontSize', 12);
ylabel(ax_ir,    'IR intensity (KM mol^{-1})',    'FontWeight', 'bold', 'FontSize', 12);
exportgraphics(fig2, 'test_2_ir_raman.pdf', 'ContentType', 'vector');
fprintf('\nDone. Figures: manual Raman plot, built-in IR+Raman plot.\n');
fprintf('Saved: test_2_raman_manual.pdf, test_2_ir_raman.pdf\n');
