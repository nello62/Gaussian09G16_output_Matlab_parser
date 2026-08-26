%% Example 4
%% Recover graphics handles from G16_draw_molecule and restyle them
close all;
clear all;
mol = G16_structure('test_2.out');
G16_draw_molecule(mol, 'Title', 'Violacein (Default Figure)');
axd  = gca;
figd = gcf;
% Default figure from G16_Draw_molecule 
exportgraphics(figd, 'fig5a.pdf', 'ContentType', 'vector');


G16_draw_molecule(mol, 'Title', 'Violacein (Default Figure)');

ax  = gca;
fig = gcf;

%% 1) Title -- always visible (HandleVisibility default 'on')
title(ax, 'Violacein (Modified Figure)', 'FontSize', 13, ...
      'FontWeight', 'bold', 'Interpreter', 'tex');

%% 2) Axis labels -- G16_draw_molecule calls axis(ax,'off'), so labels
%    exist but are not shown until the axis is switched back on
axis(ax, 'on');
xlabel(ax, 'X (\AA)', 'FontSize', 11, 'Interpreter', 'latex');
ylabel(ax, 'Y (\AA)', 'FontSize', 11, 'Interpreter', 'latex');
zlabel(ax, 'Z (\AA)', 'FontSize', 11, 'Interpreter', 'latex');

%% 3) Legend -- created with default HandleVisibility, findobj works fine
hleg = findobj(fig, 'Type', 'Legend');
set(hleg, 'FontSize', 10, 'Location', 'northeast', 'Box', 'on');

%% 4) Atom spheres (Surface objects) -- use findall: most have
%    HandleVisibility 'off' (only the first atom per element is 'on')
hatoms = findall(ax, 'Type', 'Surface');
set(hatoms, 'FaceAlpha', 0.9, 'EdgeColor', 'none');

%% 5) Bond lines -- all created with HandleVisibility 'off', findall required
hbonds = findall(ax, 'Type', 'Line');
set(hbonds, 'LineWidth', 2.5, 'Color', [0.35 0.35 0.35]);

%% 6) Atom-label text objects only, excluding Title/XLabel/YLabel/ZLabel
%    (these are also Type='text', so must be explicitly subtracted out)
htitle_h = get(ax, 'Title');
hxlab    = get(ax, 'XLabel');
hylab    = get(ax, 'YLabel');
hzlab    = get(ax, 'ZLabel');
halltext = findall(ax, 'Type', 'Text');
hatomlbl = setdiff(halltext, [htitle_h; hxlab; hylab; hzlab]);
set(hatomlbl, 'FontSize', 8, 'FontWeight', 'bold');

%% 7) Background / figure colour
set(ax,  'Color', [1 1 1]);
set(fig, 'Color', [1 1 1]);

exportgraphics(fig, 'fig5b.pdf', 'ContentType', 'vector');