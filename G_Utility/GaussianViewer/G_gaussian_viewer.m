function G_gaussian_viewer(filename)
% G_GAUSSIAN_VIEWER  GaussView-style GUI for the G09/G16 Toolbox.
%
%   G_gaussian_viewer()
%   G_gaussian_viewer(filename)
%
%   Single-window GUI wrapping the core toolbox and MOSurface: load a
%   Gaussian output file, inspect the optimised structure, browse its
%   vibrational normal modes, and render real-space MO/electron-density/
%   ESP isosurfaces (from a companion .fchk file) -- all in the same
%   embedded 3D view, so the camera orientation/zoom persist across
%   switches, in the same spirit as GaussView's own single-pane layout.
%
%   Works with BOTH Gaussian 09 and Gaussian 16 output files for the
%   structure/normal-mode part: the Gaussian major version is
%   auto-detected via G16_GAUSSIAN_VERSION (as already done by
%   G_RAMAN_BROWSER), dispatching to the matching G09_*/G16_* functions.
%   MOSurface itself (G_DRAW_MO_SURFACE/G_DRAW_DENSITY_SURFACE/
%   G_DRAW_ESP_SURFACE) is Gaussian-version-agnostic (it reads raw
%   basis-set sections directly from the .fchk file), so no dispatch is
%   needed for that tab.
%
%   Input:
%       filename - (optional) a Gaussian output file (.log/.out) to load
%                  immediately. If omitted, the viewer opens empty with
%                  an "Open file..." button.
%
%   v1 scope: structure, normal modes, and MOSurface (MO/density/ESP
%   isosurfaces). NOT included in this version (candidates for a later
%   tab, once this architecture is validated in use):
%     - IR/Raman spectrum plotting (G16_SPECTRA) and interactive
%       click-to-select Raman browsing (already available standalone as
%       G_RAMAN_BROWSER).
%     - Orbital energy-level diagram (G16_DRAW_ORBITAL).
%     - TD-DFT excited-state table, full G16_WRITE_REPORT text dump,
%       batch processing (G16_BATCH_READ_ALL).
%   ESP rendering has no background/cancel support in v1 -- it blocks the
%   UI with an indeterminate progress dialog for however long
%   G_DRAW_ESP_SURFACE itself takes (seconds to hours, depending on basis
%   set and grid resolution); a confirmation dialog warns about this
%   before starting.
%
%   Example:
%       G_gaussian_viewer();
%       G_gaussian_viewer('violacein.out');
%
%   See also G16_MODEVIEWER, G_RAMAN_BROWSER, G_DRAW_MO_SURFACE,
%            G_DRAW_DENSITY_SURFACE, G_DRAW_ESP_SURFACE.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

if nargin < 1
    filename = '';
end
if ~(isempty(filename) || ischar(filename) || (isstring(filename) && isscalar(filename)))
    error('G_gaussian_viewer:invalidInput', ...
        'G_gaussian_viewer expects a Gaussian output FILENAME (char/string), or no argument at all.');
end
filename = char(filename);

% -------------------------------------------------------------------------
% Session state (nested-function closures share these, same pattern as
% G16_MODEVIEWER / G_RAMAN_BROWSER -- no classdef app for this tool).
% -------------------------------------------------------------------------
currentFilename = '';
verLabel        = '';
mol = [];  nm = [];  en = [];  oe = [];
chg = NaN; mult = NaN;
selectedModeIdx = [];

structureFcn    = [];
nmodesFcn       = [];
drawModeFcn     = [];
drawMoleculeFcn = [];
animateModeFcn  = [];
chargeMultFcn   = [];
energyFcn       = [];
orbEnergiesFcn  = [];

fchkFile = '';
fchkData = [];

% -------------------------------------------------------------------------
% Find the active monitor, then build the window in a single atomic call.
% Same "short-lived invisible ordinary figure" trick as G16_MODEVIEWER/
% G_RAMAN_BROWSER -- reading/writing a uifigure's own Position right after
% creation is unsafe (asynchronous native-renderer handshake), so an
% ordinary FIGURE is used purely to ask MATLAB which monitor a new window
% would normally appear on.
% -------------------------------------------------------------------------
winW = 1200;
winH = 760;
tmpFig = figure('Visible', 'off');
drawnow;
tmpPos = tmpFig.Position;
delete(tmpFig);

mp = get(groot, 'MonitorPositions');
monIdx = find(tmpPos(1) >= mp(:,1) & tmpPos(1) <= mp(:,1) + mp(:,3) & ...
              tmpPos(2) >= mp(:,2) & tmpPos(2) <= mp(:,2) + mp(:,4), 1);
if isempty(monIdx), monIdx = 1; end
scr = mp(monIdx, :);

winX = scr(1) + 20;
winY = max(scr(2) + 40, scr(2) + scr(4) - winH - 80);

fig = uifigure('Name', 'G_Gaussian_Viewer', 'Position', [winX winY winW winH]);
fig.CloseRequestFcn = @(src, evt) delete(fig);

% -------------------------------------------------------------------------
% Window layout: bottom status/export strip, left sidebar, central 3D view.
% -------------------------------------------------------------------------
sidebarW = 380;

statusLabel = uilabel(fig, 'Position', [10 8 winW-410 28], ...
    'Text', 'No file loaded.', 'FontColor', [0.35 0.35 0.35]);
saveViewBtn = uibutton(fig, 'push', 'Position', [winW-330 6 150 30], ...
    'Text', 'Save view...', 'ButtonPushedFcn', @(s,e) saveView());
openBtnTop = uibutton(fig, 'push', 'Position', [winW-170 6 150 30], ...
    'Text', 'Open file...', 'ButtonPushedFcn', @(s,e) onOpenFile());

sidebar = uipanel(fig, 'Position', [0 40 sidebarW winH-40], 'BorderType', 'line');
ax = uiaxes(fig, 'Position', [sidebarW+10 45 winW-sidebarW-20 winH-50]);
ax.Toolbar.Visible = 'on';
% Interactive rotation is enabled by G09/G16_draw_molecule itself (it
% sets up rotate3d, with its own camlight-relinking callback, every time
% it renders into AX) -- no separate setup needed here before a file is
% loaded, since there is nothing to rotate yet.

% -------------------------------------------------------------------------
% Sidebar: file/summary block + tab group (Modes / MOSurface / Display)
% -------------------------------------------------------------------------
uibutton(sidebar, 'push', 'Position', [10 680 sidebarW-20 30], ...
    'Text', 'Open file...', 'FontWeight', 'bold', ...
    'ButtonPushedFcn', @(s,e) onOpenFile());

lblFile       = uilabel(sidebar, 'Position', [10 656 sidebarW-20 18], 'Text', 'File: -');
lblVersion    = uilabel(sidebar, 'Position', [10 638 sidebarW-20 18], 'Text', 'Version: -');
lblFormula    = uilabel(sidebar, 'Position', [10 620 sidebarW-20 18], 'Text', 'Formula: -');
lblChargeMult = uilabel(sidebar, 'Position', [10 602 sidebarW-20 18], 'Text', 'Charge/Mult: -');
lblEnergy     = uilabel(sidebar, 'Position', [10 584 sidebarW-20 18], 'Text', 'SCF Energy: -');
lblGap        = uilabel(sidebar, 'Position', [10 566 sidebarW-20 18], 'Text', 'HOMO-LUMO gap: -');

tg = uitabgroup(sidebar, 'Position', [5 10 sidebarW-10 540]);
tabModes   = uitab(tg, 'Title', 'Modes');
tabMO      = uitab(tg, 'Title', 'MOSurface');
tabDisplay = uitab(tg, 'Title', 'Display');

% ---- Modes tab ----------------------------------------------------------
modesTable = uitable(tabModes, 'Position', [5 45 sidebarW-30 470], ...
    'ColumnName', {'Mode','Freq (cm^-1)','Sym','IR','Raman'}, ...
    'ColumnEditable', false(1,5), ...
    'CellSelectionCallback', @(s,e) onModeSelected(e));
uibutton(tabModes, 'push', 'Position', [5 8 sidebarW-30 32], ...
    'Text', 'Animate mode (MP4)...', 'ButtonPushedFcn', @(s,e) onAnimate());

% ---- MOSurface tab -------------------------------------------------------
fchkLabel = uilabel(tabMO, 'Position', [5 495 sidebarW-30 34], ...
    'Text', 'No .fchk loaded.', 'FontColor', [0.55 0.15 0.15], 'WordWrap', 'on');
uibutton(tabMO, 'push', 'Position', [5 465 sidebarW-30 26], ...
    'Text', 'Browse for .fchk...', 'ButtonPushedFcn', @(s,e) onBrowseFchk());

uilabel(tabMO, 'Position', [5 436 sidebarW-30 18], 'Text', 'Molecular Orbital', 'FontWeight', 'bold');
uilabel(tabMO, 'Position', [5 414 200 16], 'Text', 'MO index (HOMO/LUMO/HOMO-n/n):');
moIdxField = uieditfield(tabMO, 'text', 'Position', [5 393 sidebarW-30 24], 'Value', 'HOMO');

uilabel(tabMO, 'Position', [5 364 110 18], 'Text', 'IsoValue:');
isoValField = uieditfield(tabMO, 'numeric', 'Position', [120 362 sidebarW-150 22], 'Value', 0.02);
uilabel(tabMO, 'Position', [5 336 110 18], 'Text', 'GridSpacing (Bohr):');
gridSpField = uieditfield(tabMO, 'numeric', 'Position', [120 334 sidebarW-150 22], 'Value', 0.15);
uilabel(tabMO, 'Position', [5 308 110 18], 'Text', 'Padding (Bohr):');
paddingField = uieditfield(tabMO, 'numeric', 'Position', [120 306 sidebarW-150 22], 'Value', 4.0);
showMolCheck = uicheckbox(tabMO, 'Position', [5 280 sidebarW-30 22], ...
    'Text', 'Show molecule overlay', 'Value', true);
uibutton(tabMO, 'push', 'Position', [5 248 sidebarW-30 28], ...
    'Text', 'Render MO', 'FontWeight', 'bold', 'ButtonPushedFcn', @(s,e) onRenderMO());

uilabel(tabMO, 'Position', [5 214 sidebarW-30 18], 'Text', 'Electron Density', 'FontWeight', 'bold');
spinDensityCheck = uicheckbox(tabMO, 'Position', [5 188 sidebarW-30 22], ...
    'Text', 'Spin density (UHF/UKS only)', 'Value', false);
uibutton(tabMO, 'push', 'Position', [5 156 sidebarW-30 28], ...
    'Text', 'Render Density', 'FontWeight', 'bold', 'ButtonPushedFcn', @(s,e) onRenderDensity());

uilabel(tabMO, 'Position', [5 122 sidebarW-30 18], 'Text', 'Electrostatic Potential', 'FontWeight', 'bold');
uilabel(tabMO, 'Position', [5 74 sidebarW-30 44], 'WordWrap', 'on', 'FontColor', [0.55 0.35 0.05], ...
    'Text', 'Can be slow -- from tens of seconds to several hours, depending on basis set and grid resolution. See the command window for progress.');
uibutton(tabMO, 'push', 'Position', [5 40 sidebarW-30 28], ...
    'Text', 'Render ESP...', 'FontWeight', 'bold', 'ButtonPushedFcn', @(s,e) onRenderESP());

% ---- Display tab ----------------------------------------------------------
showLabelsCheck = uicheckbox(tabDisplay, 'Position', [5 480 sidebarW-30 22], ...
    'Text', 'Show atom labels', 'Value', false);
singleBondsCheck = uicheckbox(tabDisplay, 'Position', [5 452 sidebarW-30 22], ...
    'Text', 'Single bonds only', 'Value', false);
uilabel(tabDisplay, 'Position', [5 420 150 18], 'Text', 'Bond tolerance:');
bondTolField = uieditfield(tabDisplay, 'numeric', 'Position', [160 418 sidebarW-190 22], 'Value', 1.30);
uibutton(tabDisplay, 'push', 'Position', [5 380 sidebarW-30 30], ...
    'Text', 'Redraw molecule', 'ButtonPushedFcn', @(s,e) onRedrawMolecule());

% -------------------------------------------------------------------------
if ~isempty(filename)
    loadFile(filename);
end

% =========================================================================
%  Nested callback/helper functions (share the outer function's workspace)
% =========================================================================

    function clearAx()
    % Removes every plotted object from the shared AX before a new
    % render. CLA(ax) alone does NOT fully clear a uiaxes -- unlike an
    % ordinary axes, some of its children (observed: most CPK spheres,
    % all bond lines) survive a plain CLA call, so switching between
    % molecule/mode/MOSurface renders (or between two different modes)
    % left old geometry visibly accumulating underneath the new render.
    % FINDALL(ax) recurses into every descendant regardless of nesting
    % depth (unlike ax.Children, which only lists DIRECT children and
    % was found to miss the same objects CLA does); AX itself is
    % filtered out of the result before deleting, since FINDALL(ax)
    % also returns ax.
        kids = findall(ax);
        kids(kids == ax) = [];
        delete(kids);
    end

% -------------------------------------------------------------------------
    function loadFile(fn)
    % Reads structure/modes/energy/charge/orbital-energy data from FN,
    % dispatching to G09_* or G16_* functions per the detected Gaussian
    % major version. Each secondary quantity (nm/en/oe/chg/mult) is
    % wrapped in its own try/catch so a file missing one job type (e.g.
    % no freq, no pop analysis) still loads everything else it has,
    % mirroring G16_READ_ALL's own graceful-degradation behaviour.
        if ~isfile(fn)
            uialert(fig, sprintf('File not found: %s', fn), 'Open file');
            return
        end

        try
            gv = G16_gaussian_version(fn);
        catch
            gv = struct('major', []);
        end
        if isequal(gv.major, 9)
            structureFcn    = @G09_structure;
            nmodesFcn       = @G09_nmodes;
            drawModeFcn     = @G09_draw_mode;
            drawMoleculeFcn = @G09_draw_molecule;
            animateModeFcn  = @G09_animate_mode;
            chargeMultFcn   = @G09_charge_mult;
            energyFcn       = @G09_energy;
            orbEnergiesFcn  = @G09_orbital_energies;
            verLabel        = 'Gaussian 09';
        else
            structureFcn    = @G16_structure;
            nmodesFcn       = @G16_nmodes;
            drawModeFcn     = @G16_draw_mode;
            drawMoleculeFcn = @G16_draw_molecule;
            animateModeFcn  = @G16_animate_mode;
            chargeMultFcn   = @G16_charge_mult;
            energyFcn       = @G16_energy;
            orbEnergiesFcn  = @G16_orbital_energies;
            verLabel        = 'Gaussian 16';
            if isempty(gv.major)
                verLabel = 'Gaussian 16 (assumed)';
            end
        end

        try
            mol = structureFcn(fn);
        catch ME
            uialert(fig, ME.message, 'Failed to read structure');
            return
        end

        try
            nm = nmodesFcn(fn);
        catch
            nm = [];
        end
        try
            [chg, mult] = chargeMultFcn(fn);
        catch
            chg = NaN; mult = NaN;
        end
        try
            en = energyFcn(fn);
        catch
            en = [];
        end
        try
            oe = orbEnergiesFcn(fn);
        catch
            oe = [];
        end

        currentFilename = fn;
        selectedModeIdx = [];

        updateSummary();
        updateModesTable();
        redrawMolecule();

        % Auto-detect a sibling .fchk (same basename) for the MOSurface tab.
        [p_, b_, ~] = fileparts(fn);
        candidate = fullfile(p_, [b_, '.fchk']);
        if isfile(candidate)
            setFchk(candidate);
        else
            fchkData = [];
            fchkFile = '';
            fchkLabel.Text = 'No .fchk loaded (browse to enable MOSurface).';
            fchkLabel.FontColor = [0.55 0.15 0.15];
        end
    end

% -------------------------------------------------------------------------
    function onOpenFile()
        [f, p] = uigetfile({'*.out;*.log', 'Gaussian output (*.out, *.log)'; ...
                             '*.fchk', 'Formatted checkpoint (*.fchk)'}, ...
                             'Open Gaussian file');
        if isequal(f, 0)
            return
        end
        loadFile(fullfile(p, f));
    end

% -------------------------------------------------------------------------
    function updateSummary()
        [~, fname_] = fileparts(currentFilename);
        lblFile.Text    = sprintf('File: %s', fname_);
        lblVersion.Text = sprintf('Version: %s', verLabel);
        lblFormula.Text = sprintf('Formula: %s (%d atoms)', formula_from_symbols(mol.symbols), mol.Natoms);
        if isnan(chg)
            lblChargeMult.Text = 'Charge/Mult: -';
        else
            lblChargeMult.Text = sprintf('Charge/Mult: %+d / %d', chg, mult);
        end
        if ~isempty(en)
            lblEnergy.Text = sprintf('SCF Energy: %.6f Ha (%s)', en.SCF, en.method);
        else
            lblEnergy.Text = 'SCF Energy: n/a';
        end
        if ~isempty(oe)
            lblGap.Text = sprintf('HOMO-LUMO gap: %.3f eV', oe.gap_eV);
        else
            lblGap.Text = 'HOMO-LUMO gap: n/a';
        end
    end

% -------------------------------------------------------------------------
    function updateModesTable()
        if isempty(nm) || nm.Nmodes < 1
            modesTable.Data = {};
            return
        end
        data = cell(nm.Nmodes, 5);
        for i = 1:nm.Nmodes
            symLabel = '';
            if ~isempty(nm.symmetry) && numel(nm.symmetry) >= i && ~isempty(nm.symmetry{i})
                symLabel = nm.symmetry{i};
            end
            data{i,1} = i;
            data{i,2} = nm.freq(i);
            data{i,3} = symLabel;
            data{i,4} = nm.IR(i);
            if nm.has_Raman
                data{i,5} = nm.Raman(i);
            else
                data{i,5} = NaN;
            end
        end
        modesTable.Data = data;
    end

% -------------------------------------------------------------------------
    function redrawMolecule()
        if isempty(mol)
            return
        end
        clearAx();
        try
            drawMoleculeFcn(mol, 'Ax', ax, 'ShowLabels', showLabelsCheck.Value, ...
                'SingleBondsOnly', singleBondsCheck.Value, 'BondTol', bondTolField.Value, ...
                'Title', '');
            [~, fname_] = fileparts(currentFilename);
            statusLabel.Text = sprintf('Molecule: %s', fname_);
        catch ME
            uialert(fig, ME.message, 'draw_molecule error');
        end
    end

% -------------------------------------------------------------------------
    function onRedrawMolecule()
        redrawMolecule();
    end

% -------------------------------------------------------------------------
    function onModeSelected(evt)
        if isempty(evt.Indices)
            return
        end
        row = evt.Indices(1,1);
        idx = modesTable.Data{row,1};
        selectedModeIdx = idx;
        clearAx();
        try
            drawModeFcn(mol, nm, idx, 'Ax', ax, 'ShowLabels', showLabelsCheck.Value, ...
                'BondTol', bondTolField.Value);
            statusLabel.Text = sprintf('Mode %d -- %.1f cm^{-1}', idx, nm.freq(idx));
        catch ME
            uialert(fig, ME.message, 'draw_mode error');
        end
    end

% -------------------------------------------------------------------------
    function onAnimate()
        if isempty(selectedModeIdx)
            uialert(fig, 'Select a mode from the table first.', 'No mode selected');
            return
        end
        [~, base_] = fileparts(currentFilename);
        defaultName = sprintf('%s_mode%d.mp4', matlab.lang.makeValidName(base_), selectedModeIdx);
        [f, p] = uiputfile({'*.mp4', 'MP4 video'}, 'Save mode animation as', defaultName);
        if isequal(f, 0)
            return
        end
        outFile = fullfile(p, f);

        [az, el] = view(ax);
        d = uiprogressdlg(fig, 'Title', 'Rendering animation...', 'Indeterminate', 'on');
        try
            animateModeFcn(mol, nm, selectedModeIdx, 'Filename', outFile, 'View', [az el]);
            statusLabel.Text = sprintf('Animation saved: %s', outFile);
        catch ME
            uialert(fig, ME.message, 'animate_mode error');
        end
        close(d);
    end

% -------------------------------------------------------------------------
    function setFchk(fn)
        try
            fchkData = G16_fchk_read(fn);
            fchkFile = fn;
            [~, b_] = fileparts(fn);
            fchkLabel.Text = sprintf('.fchk: %s', b_);
            fchkLabel.FontColor = [0.10 0.45 0.10];
        catch ME
            fchkData = [];
            fchkFile = '';
            fchkLabel.Text = 'Failed to read .fchk (see error).';
            fchkLabel.FontColor = [0.55 0.15 0.15];
            uialert(fig, ME.message, 'G16_fchk_read error');
        end
    end

% -------------------------------------------------------------------------
    function onBrowseFchk()
        [f, p] = uigetfile({'*.fchk', 'Formatted checkpoint (*.fchk)'}, 'Select .fchk file');
        if isequal(f, 0)
            return
        end
        setFchk(fullfile(p, f));
    end

% -------------------------------------------------------------------------
    function moArg = resolveMoField()
    % MO index field accepts 'HOMO'/'LUMO'/'HOMO-n'/'LUMO+n' (forwarded
    % as a string) or a plain integer (forwarded as a number) --
    % G_draw_mo_surface's own resolve_mo_index validates either form, no
    % need to duplicate that validation here.
        s = strtrim(moIdxField.Value);
        v = str2double(s);
        if isnan(v)
            moArg = s;
        else
            moArg = v;
        end
    end

% -------------------------------------------------------------------------
    function onRenderMO()
        if isempty(fchkData)
            uialert(fig, 'Load a .fchk file first (MOSurface tab).', 'No .fchk loaded');
            return
        end
        clearAx();
        d = uiprogressdlg(fig, 'Title', 'Rendering molecular orbital...', 'Indeterminate', 'on');
        try
            G_draw_mo_surface(fchkData, resolveMoField(), 'Ax', ax, ...
                'IsoValue', isoValField.Value, 'GridSpacing', gridSpField.Value, ...
                'Padding', paddingField.Value, 'ShowMolecule', showMolCheck.Value);
            statusLabel.Text = sprintf('MO surface: %s', moIdxField.Value);
        catch ME
            uialert(fig, ME.message, 'G_draw_mo_surface error');
        end
        close(d);
    end

% -------------------------------------------------------------------------
    function onRenderDensity()
        if isempty(fchkData)
            uialert(fig, 'Load a .fchk file first (MOSurface tab).', 'No .fchk loaded');
            return
        end
        clearAx();
        d = uiprogressdlg(fig, 'Title', 'Rendering electron density...', 'Indeterminate', 'on');
        try
            G_draw_density_surface(fchkData, 'Ax', ax, 'GridSpacing', gridSpField.Value, ...
                'Padding', paddingField.Value, 'ShowMolecule', showMolCheck.Value, ...
                'SpinDensity', spinDensityCheck.Value);
            if spinDensityCheck.Value
                statusLabel.Text = 'Spin density surface';
            else
                statusLabel.Text = 'Electron density surface';
            end
        catch ME
            uialert(fig, ME.message, 'G_draw_density_surface error');
        end
        close(d);
    end

% -------------------------------------------------------------------------
    function onRenderESP()
        if isempty(fchkData)
            uialert(fig, 'Load a .fchk file first (MOSurface tab).', 'No .fchk loaded');
            return
        end
        sel = uiconfirm(fig, ...
            ['ESP evaluation can take anywhere from tens of seconds to several ' ...
             'hours, depending on the basis set and grid resolution -- see the ' ...
             'command window for progress once it starts. Continue?'], ...
            'Render ESP surface', 'Options', {'Run', 'Cancel'}, ...
            'DefaultOption', 2, 'CancelOption', 2);
        if ~strcmp(sel, 'Run')
            return
        end
        clearAx();
        d = uiprogressdlg(fig, 'Title', 'Rendering ESP surface...', ...
            'Message', 'This may take a while -- see the command window for progress.', ...
            'Indeterminate', 'on');
        try
            G_draw_esp_surface(fchkData, 'Ax', ax, 'GridSpacing', gridSpField.Value, ...
                'Padding', paddingField.Value, 'ShowMolecule', showMolCheck.Value);
            statusLabel.Text = 'ESP surface';
        catch ME
            uialert(fig, ME.message, 'G_draw_esp_surface error');
        end
        close(d);
    end

% -------------------------------------------------------------------------
    function saveView()
        if isempty(mol)
            uialert(fig, 'Nothing to save -- load a file first.', 'Save view');
            return
        end
        [~, base_] = fileparts(currentFilename);
        defaultName = [matlab.lang.makeValidName(base_), '.pdf'];
        [f, p] = uiputfile({'*.pdf', 'PDF (vector)'; '*.eps', 'EPS (vector)'; ...
                             '*.png', 'PNG (raster)'; '*.jpg', 'JPEG (raster)'}, ...
                             'Save current view as', defaultName);
        if isequal(f, 0)
            return
        end
        outFile = fullfile(p, f);
        [~, ~, ext] = fileparts(outFile);
        try
            if any(strcmpi(ext, {'.pdf', '.eps'}))
                exportgraphics(ax, outFile, 'ContentType', 'vector');
            else
                exportgraphics(ax, outFile, 'Resolution', 300);
            end
            statusLabel.Text = sprintf('View saved: %s', outFile);
        catch ME
            uialert(fig, ME.message, 'Export error');
        end
    end

end % function G_gaussian_viewer


% =========================================================================
%  Local function: element-count formula string (Hill system: C, H, then
%  the remaining elements alphabetically; or fully alphabetical if no C).
% =========================================================================
function s = formula_from_symbols(symbols)
counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:numel(symbols)
    sym = symbols{i};
    if isKey(counts, sym)
        counts(sym) = counts(sym) + 1;
    else
        counts(sym) = 1;
    end
end
elems = keys(counts);
if any(strcmp(elems, 'C'))
    order = {'C'};
    if any(strcmp(elems, 'H'))
        order{end+1} = 'H';
    end
    order = [order, sort(setdiff(elems, {'C','H'}))];
else
    order = sort(elems);
end
parts = cell(1, numel(order));
for k = 1:numel(order)
    n = counts(order{k});
    if n == 1
        parts{k} = order{k};
    else
        parts{k} = sprintf('%s%d', order{k}, n);
    end
end
s = strjoin(parts, '');
end
