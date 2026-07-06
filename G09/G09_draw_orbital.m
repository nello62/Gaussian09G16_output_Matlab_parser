function G09_draw_orbital(oe, varargin)
% G09_DRAW_ORBITAL  Draws a molecular-orbital energy-level diagram from
%                    the struct returned by G09_ORBITAL_ENERGIES,
%                    highlighting the HOMO-LUMO transition with an arrow
%                    labelled with the gap.
%
%   G09_draw_orbital(oe)
%   G09_draw_orbital(oe, Name, Value, ...)
%
%   Input:
%       oe  - struct returned by G09_orbital_energies
%
%   Optional parameters (Name-Value):
%       'NLevels'     - number of occupied/virtual levels to show around
%                       the frontier orbitals, on each side (default: 5)
%       'Units'       - 'eV' (default) | 'Hartree'
%       'Ax'          - existing axes handle (default: new figure)
%       'OccColor'    - colour for non-frontier occupied levels (default: [0.25 0.25 0.70])
%       'VirtColor'   - colour for non-frontier virtual levels (default: [0.70 0.30 0.25])
%       'HOMOColor'   - highlight colour for HOMO (default: [0 0.45 0.85])
%       'LUMOColor'   - highlight colour for LUMO (default: [0.90 0.35 0])
%       'ArrowColor'  - HOMO-LUMO transition arrow colour (default: [0.15 0.60 0.15])
%       'ShowLabels'  - annotate each level with its energy value (default: true)
%       'ShowSpins'   - draw a paired-electron arrow glyph on occupied levels (default: true)
%       'FontSize'    - level energy-value label font size (default: 8)
%       'FrontierFontSize' - font size for the HOMO/LUMO tags and the
%                       gap annotation (default: 13)
%       'ArrowHeadSize' - HOMO-LUMO arrow head size, as a fraction of the
%                       arrow length, passed to QUIVER (default: 0.2)
%       'Title'       - figure title (default: filename or 'Orbital Energy Diagram')
%
%   Example:
%       oe = G09_orbital_energies('a1.out');
%       G09_draw_orbital(oe)
%       G09_draw_orbital(oe, 'NLevels', 8, 'Units', 'Hartree')
%
%   See also G09_ORBITAL_ENERGIES.

% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'oe');
addParameter(p, 'NLevels',    5,                @isnumeric);
addParameter(p, 'Units',      'eV',             @ischar);
addParameter(p, 'Ax',         [],               @ishandle);
addParameter(p, 'OccColor',   [0.25 0.25 0.70], @isnumeric);
addParameter(p, 'VirtColor',  [0.70 0.30 0.25], @isnumeric);
addParameter(p, 'HOMOColor',  [0.00 0.45 0.85], @isnumeric);
addParameter(p, 'LUMOColor',  [0.90 0.35 0.00], @isnumeric);
addParameter(p, 'ArrowColor', [0.15 0.60 0.15], @isnumeric);
addParameter(p, 'ShowLabels', true,             @islogical);
addParameter(p, 'ShowSpins',  true,             @islogical);
addParameter(p, 'FontSize',   8,                @isnumeric);
addParameter(p, 'FrontierFontSize', 13,         @isnumeric);
addParameter(p, 'ArrowHeadSize', 0.2,           @isnumeric);
addParameter(p, 'Title',      '',               @ischar);
parse(p, oe, varargin{:});

nlevels     = round(p.Results.NLevels);
units       = lower(p.Results.Units);
ax          = p.Results.Ax;
occ_color   = p.Results.OccColor;
virt_color  = p.Results.VirtColor;
homo_color  = p.Results.HOMOColor;
lumo_color  = p.Results.LUMOColor;
arrow_color = p.Results.ArrowColor;
show_labels = p.Results.ShowLabels;
show_spins  = p.Results.ShowSpins;
fsize       = p.Results.FontSize;
fsize_front = p.Results.FrontierFontSize;
arrow_head  = p.Results.ArrowHeadSize;
fig_title   = p.Results.Title;

% Validate oe struct
if ~isstruct(oe) || ~isfield(oe, 'alpha_occ') || ~isfield(oe, 'alpha_virt')
    error('G09_draw_orbital: oe must be the struct returned by G09_orbital_energies.');
end
if isempty(oe.alpha_occ) || isempty(oe.alpha_virt)
    error('G09_draw_orbital: oe has no occupied/virtual orbital energies to plot.');
end

% -------------------------------------------------------------------------
% Unit conversion
% -------------------------------------------------------------------------
ha2eV = 27.211386245988;
switch units
    case 'ev'
        occ_e    = oe.alpha_occ  * ha2eV;
        virt_e   = oe.alpha_virt * ha2eV;
        gap_disp = oe.gap_eV;
        unit_lbl = 'eV';
    case 'hartree'
        occ_e    = oe.alpha_occ;
        virt_e   = oe.alpha_virt;
        gap_disp = oe.gap;
        unit_lbl = 'Ha';
    otherwise
        error('G09_draw_orbital: ''Units'' must be ''eV'' or ''Hartree''.');
end

Nocc  = numel(occ_e);
Nvirt = numel(virt_e);
occ_show  = occ_e(max(1, Nocc-nlevels+1):Nocc);   % last N occupied (HOMO last), ascending energy
virt_show = virt_e(1:min(nlevels, Nvirt));         % first N virtual (LUMO first), ascending energy

HOMO_e = occ_e(end);
LUMO_e = virt_e(1);

% -------------------------------------------------------------------------
% Figure and axes
% -------------------------------------------------------------------------
if isempty(ax)
    fig = figure('Color', 'white', 'Name', 'Orbital energy diagram', 'NumberTitle', 'off');
    ax  = axes('Parent', fig);
end
hold(ax, 'on');

bar_x = [-0.5, 0.5];

% -------------------------------------------------------------------------
% Occupied levels
% -------------------------------------------------------------------------
for i = 1:numel(occ_show)
    e = occ_show(i);
    is_homo = (i == numel(occ_show));
    if is_homo
        clr = homo_color; lw = 2.5;
    else
        clr = occ_color; lw = 1.5;
    end
    line(ax, bar_x, [e e], 'Color', clr, 'LineWidth', lw, 'HandleVisibility', 'off');
    if show_spins
        text(ax, 0, e, '\uparrow\downarrow', 'FontSize', fsize+2, 'Color', clr, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'Interpreter', 'tex', 'HandleVisibility', 'off');
    end
    if show_labels
        text(ax, bar_x(1)-0.08, e, sprintf('%.3f', e), 'FontSize', fsize, 'Color', clr, ...
             'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
             'Interpreter', 'none', 'HandleVisibility', 'off');
    end
end

% -------------------------------------------------------------------------
% Virtual levels
% -------------------------------------------------------------------------
for i = 1:numel(virt_show)
    e = virt_show(i);
    is_lumo = (i == 1);
    if is_lumo
        clr = lumo_color; lw = 2.5;
    else
        clr = virt_color; lw = 1.5;
    end
    line(ax, bar_x, [e e], 'Color', clr, 'LineWidth', lw, 'HandleVisibility', 'off');
    if show_labels
        text(ax, bar_x(1)-0.08, e, sprintf('%.3f', e), 'FontSize', fsize, 'Color', clr, ...
             'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
             'Interpreter', 'none', 'HandleVisibility', 'off');
    end
end

% -------------------------------------------------------------------------
% HOMO/LUMO tags and HOMO -> LUMO transition arrow, with the gap value
% -------------------------------------------------------------------------
text(ax, bar_x(2)+0.08, HOMO_e, 'HOMO', 'Color', homo_color, 'FontWeight', 'bold', ...
     'FontSize', fsize_front, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
     'Interpreter', 'none', 'HandleVisibility', 'off');
text(ax, bar_x(2)+0.08, LUMO_e, 'LUMO', 'Color', lumo_color, 'FontWeight', 'bold', ...
     'FontSize', fsize_front, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
     'Interpreter', 'none', 'HandleVisibility', 'off');

arrow_x = 0.95;
quiver(ax, arrow_x, HOMO_e, 0, LUMO_e-HOMO_e, 0, 'Color', arrow_color, ...
    'LineWidth', 2, 'MaxHeadSize', arrow_head, 'HandleVisibility', 'off');
text(ax, arrow_x+0.12, (HOMO_e+LUMO_e)/2, sprintf('\\DeltaE = %.3f %s', gap_disp, unit_lbl), ...
     'FontSize', fsize_front, 'FontWeight', 'bold', 'Color', arrow_color, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
     'Interpreter', 'tex', 'HandleVisibility', 'off');

% -------------------------------------------------------------------------
% Cosmetics
% -------------------------------------------------------------------------
xlim(ax, [-1.0, 2.3]);
set(ax, 'XTick', [], 'Box', 'off');
ylabel(ax, sprintf('Energy (%s)', unit_lbl));

if isempty(fig_title)
    if isfield(oe, 'filename') && ~isempty(oe.filename)
        [~, fn]   = fileparts(oe.filename);
        fig_title = strrep(fn, '_', '\_');
    else
        fig_title = 'Orbital Energy Diagram';
    end
end
title(ax, fig_title, 'Interpreter', 'tex', 'FontSize', 11);

hold(ax, 'off');

end % G09_draw_orbital
