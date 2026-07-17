function G09_modeViewer(filename, varargin)
%G09_MODEVIEWER  Interactive combo-box selector for Gaussian vibrational modes.
%
%   G09_MODEVIEWER(filename) reads the molecular structure and all
%   vibrational normal modes from FILENAME using G09_STRUCTURE and
%   G09_NMODES, then opens a selector window listing every mode (index,
%   frequency, symmetry). Selecting an entry from the drop-down calls
%   G09_DRAW_MODE to render that mode's CPK structure with displacement
%   arrows in a NEW figure window. Previously drawn mode figures are left
%   open, so different modes (or different render settings) can be
%   compared side by side.
%
%   The window also provides:
%     - An "Order by:" combo box to sort the mode list by mode number
%       (default), IR intensity, or Raman intensity (descending; only
%       offered when Raman data is present). The currently selected mode
%       stays selected across reordering.
%     - A text box to set a custom title on the current mode figure.
%     - A combo box / edit field for every G09_DRAW_MODE option (Scale,
%       ArrowColor, AtomScale, BondTol, ShowLabels, FlipSign). Changing
%       any of them immediately redraws the currently selected mode with
%       the new setting, in a new figure.
%     - A "Save figure..." button exporting the target mode figure to
%       JPEG, EPS, or PDF (EPS/PDF as true vector graphics).
%
%   G09_MODEVIEWER(filename, 'Name', Value, ...) pre-populates the
%   corresponding option control(s) with the given defaults, e.g.:
%
%       G09_modeViewer('INDACO_E025_ppDP.LOG', 'Scale', 2, 'ShowLabels', true)
%
%   Example:
%       G09_modeViewer('V_E00t.out');
%
%   See also G09_STRUCTURE, G09_NMODES, G09_DRAW_MODE.

% -------------------------------------------------------------------------
% ---- load structure and normal modes -----------------------------------
fprintf('G09_modeViewer: reading structure and normal modes from %s ...\n', filename);
mol = G09_structure(filename);
nm  = G09_nmodes(filename);
fprintf('  %d atoms, %d vibrational modes.\n', nm.Natoms, nm.Nmodes);

if nm.Nmodes < 1
    error('G09_modeViewer:noModes', 'No vibrational modes found in %s.', filename);
end

currentModeFig = gobjects(0);   % handle of the most recently drawn mode figure
openModeFigs   = gobjects(0);   % all mode figures opened by this viewer

% ---- screen geometry, so mode figures can be opened as large as possible -
% NOTE: get(groot,'ScreenSize') returns the bounding box of ALL monitors
% combined on multi-display systems (a well-known MATLAB gotcha), which can
% place windows in a region not covered by any single physical screen.
% MonitorPositions gives one row per real monitor; row 1 is normally the
% primary display, but we explicitly pick the one that contains the origin
% (1,1) to be safe regardless of monitor order/arrangement.
mp  = get(groot, 'MonitorPositions');       % [Nmonitors x 4], each row [left bottom width height]
priIdx = find(mp(:,1) <= 1 & mp(:,2) <= 1, 1);
if isempty(priIdx), priIdx = 1; end
scr = mp(priIdx, :);                        % [left bottom width height] px, primary monitor only
panelW = 480;                              % width reserved for the selector panel
panelH = 556;
selX = scr(1) + 20;
selY = max(scr(2) + 40, scr(2) + scr(4) - panelH - 80);   % near top-left, leaving room for title/taskbars

% ---- named arrow-colour presets -----------------------------------------
colorNames = {'Orange (default)','Red','Blue','Green','Black','Magenta','Custom...'};
colorMap   = containers.Map(colorNames(1:end-1), ...
    {[1 0.4 0.1], [0.85 0 0], [0 0.35 0.85], [0.10 0.60 0.20], [0 0 0], [0.75 0 0.75]});
customColor = [1 0.4 0.1];   % updated if the user picks a custom colour

% ---- parse incoming defaults / pass-through options ----------------------
defScale      = 1.5;
defAtomScale  = 0.35;
defBondTol    = 1.30;
defShowLabels = false;
defFontSize   = 10;   % post-processing only: not a G09_draw_mode parameter
defArrowColorName = 'Orange (default)';
extraOpts = {};   % any other Name-Value pair is forwarded as-is, unmodified

k = 1;
while k <= numel(varargin)-1
    name = varargin{k}; val = varargin{k+1};
    switch lower(name)
        case 'scale',      defScale = val;
        case 'atomscale',  defAtomScale = val;
        case 'bondtol',    defBondTol = val;
        case 'showlabels', defShowLabels = logical(val);
        case 'fontsize',   defFontSize = val;
        case 'flipsign',   defFlipSign = logical(val); %#ok<NASGU>
        case 'arrowcolor'
            customColor = val;
            defArrowColorName = 'Custom...';
        otherwise
            extraOpts = [extraOpts, {name, val}]; %#ok<AGROW>
    end
    k = k + 2;
end
if ~exist('defFlipSign', 'var'), defFlipSign = false; end

% ---- build drop-down entries: "Mode  12    456.3 cm-1  (A)" -------------
% ---- ordering options for the mode drop-down (mode number / IR / Raman) -
orderNames = {'Mode number', 'IR intensity'};
if nm.has_Raman
    orderNames{end+1} = 'Raman intensity';
end
[items, itemsData] = buildModeItems('Mode number');

fmtItems     = {'PDF (vector)', 'EPS (vector)', 'JPEG (raster)'};
fmtItemsData = {'pdf', 'epsc', 'jpeg'};

% ---- build the selector window -----------------------------------------
selFig = uifigure('Name', sprintf('G09 Normal Mode Viewer - %s', filename), ...
                   'Position', [selX selY 460 556]);
selFig.CloseRequestFcn = @(src, evt) closeViewer();

uilabel(selFig, 'Position', [20 514 420 22], ...
    'Text', 'Select a vibrational mode:', 'FontWeight', 'bold');

modeDD = uidropdown(selFig, ...
    'Position', [20 481 420 28], ...
    'Items', items, 'ItemsData', itemsData, 'Value', itemsData(1), ...
    'ValueChangedFcn', @(src, evt) drawMode(src.Value, false));

uilabel(selFig, 'Position', [20 421 70 18], 'Text', 'Order by:');
orderDD = uidropdown(selFig, 'Position', [95 418 200 28], ...
    'Items', orderNames, 'Value', 'Mode number', ...
    'ValueChangedFcn', @(src, evt) onOrderChanged(src.Value));

uilabel(selFig, 'Position', [20 384 420 22], ...
    'Text', 'Figure title:', 'FontWeight', 'bold');

titleEdit = uieditfield(selFig, 'text', ...
    'Position', [20 351 300 28], ...
    'Placeholder', 'e.g. Mode 70 - C=C stretch', ...
    'ValueChangedFcn', @(src, evt) applyTitle());

uibutton(selFig, 'push', 'Position', [330 351 110 28], ...
    'Text', 'Apply title', 'ButtonPushedFcn', @(src, evt) applyTitle());

uilabel(selFig, 'Position', [20 314 420 22], ...
    'Text', 'G09_draw_mode options:', 'FontWeight', 'bold');

% -- Row 1: Scale / ArrowColor --
uilabel(selFig, 'Position', [20 281 60 18], 'Text', 'Scale:');
scaleDD = uidropdown(selFig, 'Position', [85 278 130 28], 'Editable', 'on', ...
    'Items', {'0.5','1','1.5','2','3'}, 'Value', num2str(defScale), ...
    'ValueChangedFcn', @(src, evt) onOptionChanged());

uilabel(selFig, 'Position', [240 281 85 18], 'Text', 'ArrowColor:');
colorDD = uidropdown(selFig, 'Position', [330 278 110 28], ...
    'Items', colorNames, 'Value', defArrowColorName, ...
    'ValueChangedFcn', @(src, evt) onColorChanged(src));

% -- Row 2: AtomScale / BondTol --
uilabel(selFig, 'Position', [20 245 75 18], 'Text', 'AtomScale:');
atomScaleDD = uidropdown(selFig, 'Position', [100 242 115 28], 'Editable', 'on', ...
    'Items', {'0.2','0.3','0.35','0.5'}, 'Value', num2str(defAtomScale), ...
    'ValueChangedFcn', @(src, evt) onOptionChanged());

uilabel(selFig, 'Position', [240 245 65 18], 'Text', 'BondTol:');
bondTolDD = uidropdown(selFig, 'Position', [310 242 130 28], 'Editable', 'on', ...
    'Items', {'1.10','1.20','1.30','1.40','1.50'}, 'Value', num2str(defBondTol), ...
    'ValueChangedFcn', @(src, evt) onOptionChanged());

% -- Row 3: ShowLabels / FlipSign --
uilabel(selFig, 'Position', [20 209 85 18], 'Text', 'ShowLabels:');
showLabelsDD = uidropdown(selFig, 'Position', [110 206 110 28], ...
    'Items', {'false','true'}, 'Value', tf2str(defShowLabels), ...
    'ValueChangedFcn', @(src, evt) onOptionChanged());

uilabel(selFig, 'Position', [240 209 65 18], 'Text', 'FlipSign:');
flipSignDD = uidropdown(selFig, 'Position', [310 206 130 28], ...
    'Items', {'false','true'}, 'Value', tf2str(defFlipSign), ...
    'ValueChangedFcn', @(src, evt) onOptionChanged());

% -- Row 4: FontSize of the atom-index labels (post-processing; not a
%    G09_draw_mode parameter - applied to the label text objects after
%    the figure is drawn). Only has a visible effect when ShowLabels=true.
uilabel(selFig, 'Position', [20 173 70 18], 'Text', 'LabelFontSize:');
fontSizeDD = uidropdown(selFig, 'Position', [115 170 80 28], 'Editable', 'on', ...
    'Items', {'6','8','10','12','14','16'}, 'Value', num2str(defFontSize), ...
    'ValueChangedFcn', @(src, evt) onOptionChanged());

uilabel(selFig, 'Position', [20 133 420 22], ...
    'Text', 'Save current mode figure as:', 'FontWeight', 'bold');

fmtDD = uidropdown(selFig, 'Position', [20 100 180 28], ...
    'Items', fmtItems, 'ItemsData', fmtItemsData, 'Value', 'pdf');

uibutton(selFig, 'push', 'Position', [215 100 205 28], ...
    'Text', 'Save figure...', 'ButtonPushedFcn', @(src, evt) saveCurrentFigure(fmtDD.Value));

uibutton(selFig, 'push', 'Position', [20 55 420 30], ...
    'Text', 'Animate mode (MP4)...', 'FontWeight', 'bold', ...
    'ButtonPushedFcn', @(src, evt) animateCurrentMode());

uilabel(selFig, 'Position', [20 15 420 30], ...
    'Text', sprintf('%d atoms  |  %d modes  |  file: %s', nm.Natoms, nm.Nmodes, filename), ...
    'FontColor', [0.45 0.45 0.45], 'FontSize', 11);

% draw the first mode as soon as the viewer opens
drawMode(itemsData(1), false);

% =========================================================================
    function opts = buildDrawOpts()
    % Reads all option controls and returns the Name-Value cell array to
    % forward to G09_draw_mode, plus any pass-through options given at
    % construction time.
        scaleVal     = str2double(scaleDD.Value);
        atomScaleVal = str2double(atomScaleDD.Value);
        bondTolVal   = str2double(bondTolDD.Value);
        showLabelsVal = strcmp(showLabelsDD.Value, 'true');
        flipSignVal   = strcmp(flipSignDD.Value, 'true');

        if strcmp(colorDD.Value, 'Custom...')
            arrowColorVal = customColor;
        else
            arrowColorVal = colorMap(colorDD.Value);
        end

        opts = [{'Scale', scaleVal, 'ArrowColor', arrowColorVal, ...
                 'AtomScale', atomScaleVal, 'BondTol', bondTolVal, ...
                 'ShowLabels', showLabelsVal, 'FlipSign', flipSignVal}, extraOpts];
    end

% -------------------------------------------------------------------------
    function [modeItems, modeItemsData] = buildModeItems(orderName)
    % Builds the drop-down entries and their underlying mode indices,
    % ordered by mode number (ascending) or by IR/Raman intensity
    % (descending). ItemsData always holds the true mode index, so
    % reordering never changes which mode a given entry draws.
        switch orderName
            case 'IR intensity'
                [~, idx] = sort(nm.IR, 'descend');
            case 'Raman intensity'
                [~, idx] = sort(nm.Raman, 'descend');
            otherwise   % 'Mode number'
                idx = (1:nm.Nmodes)';
        end

        modeItems = cell(nm.Nmodes, 1);
        for kk = 1:nm.Nmodes
            i = idx(kk);
            symLabel = '';
            if ~isempty(nm.symmetry) && numel(nm.symmetry) >= i && ~isempty(nm.symmetry{i})
                symLabel = sprintf('  (%s)', nm.symmetry{i});
            end
            flag = '';
            if nm.freq(i) < 0
                flag = '   [imaginary]';
            end
            modeItems{kk} = sprintf('Mode %3d   %9.1f cm^-1%s%s', i, nm.freq(i), symLabel, flag);
            switch orderName
                case 'IR intensity'
                    modeItems{kk} = sprintf('%s   IR=%.1f', modeItems{kk}, nm.IR(i));
                case 'Raman intensity'
                    modeItems{kk} = sprintf('%s   Raman=%.1f', modeItems{kk}, nm.Raman(i));
            end
        end
        modeItemsData = idx;
    end

% -------------------------------------------------------------------------
    function onOrderChanged(orderName)
    % Reorders the mode drop-down entries in place, preserving whichever
    % mode is currently selected (ItemsData still maps to the same mode
    % indices, only their order/labels change).
        [newItems, newData] = buildModeItems(orderName);
        curVal = modeDD.Value;
        set(modeDD, 'Items', newItems, 'ItemsData', newData);
        modeDD.Value = curVal;
    end

% -------------------------------------------------------------------------
    function drawMode(k, closeOld)
    % Draws mode K via G09_draw_mode.
    %   closeOld = false : used when the MODE selection changes -> opens a
    %              NEW figure and keeps all previous ones open, so
    %              different modes can be compared side by side.
    %   closeOld = true  : used when a G09_draw_mode OPTION changes
    %              (Scale, ArrowColor, AtomScale, BondTol, ShowLabels,
    %              FlipSign) -> the newly rendered figure replaces the
    %              current one (old one is closed AFTER the new one is
    %              successfully drawn, to avoid any close/redraw race).
        toClose = gobjects(0);
        if closeOld && ~isempty(currentModeFig) && isgraphics(currentModeFig)
            toClose = currentModeFig;
        end

        figsBefore = findall(0, 'Type', 'figure');

        try
            optsCell = buildDrawOpts();
            G09_draw_mode(mol, nm, k, optsCell{:});
        catch ME
            uialert(selFig, ME.message, 'G09_draw_mode error');
            return
        end

        figsAfter = findall(0, 'Type', 'figure');
        newFigs   = setdiff(figsAfter, figsBefore);

        if ~isempty(newFigs)
            currentModeFig = newFigs(1);
        else
            % G09_draw_mode drew into the already-current figure
            currentModeFig = gcf;
        end
        currentModeFig.Name = sprintf('%s - Mode %d (%.1f cm^-1)', filename, k, nm.freq(k));

        % Apply the requested label font size (post-processing: this is
        % not a G09_draw_mode parameter, it edits the label text objects
        % directly after the figure has been drawn).
        applyLabelFontSize(currentModeFig, str2double(fontSizeDD.Value));

        % Make the mode figure as large as the screen allows, positioned
        % to the right of the selector panel. When several mode figures
        % are kept open at once (mode comparison), cascade them slightly
        % so each one stays visible instead of exactly overlapping.
        cascadeStep = 40;
        cascadeMax  = 6;
        nOpenBefore = numel(openModeFigs(isgraphics(openModeFigs)));
        idx = mod(nOpenBefore, cascadeMax);

        figW = max(500, scr(3) - panelW - 40 - cascadeMax*cascadeStep);
        figH = max(400, scr(4) - 80 - cascadeMax*cascadeStep);
        figX = scr(1) + panelW + idx*cascadeStep;
        figY = scr(2) + 40 + (cascadeMax - idx)*cascadeStep;
        try
            set(currentModeFig, 'Units', 'pixels', 'Position', [figX figY figW figH]);
        catch
            % ignore if the figure type does not support Position (unlikely)
        end

        % close the previous figure only now that the new one is up
        if ~isempty(toClose) && isgraphics(toClose) && toClose ~= currentModeFig
            close(toClose);
            openModeFigs = openModeFigs(openModeFigs ~= toClose);
        end

        openModeFigs = openModeFigs(isgraphics(openModeFigs));  % drop stale handles
        openModeFigs(end+1) = currentModeFig; %#ok<AGROW>
    end

% -------------------------------------------------------------------------
    function applyLabelFontSize(fig, fsz)
    % Sets the FontSize of the atom-index label text objects drawn by
    % G09_draw_mode (only present when ShowLabels = true). G09_draw_mode
    % has no FontSize parameter of its own, so this edits the text
    % objects directly, in place, after the figure has been rendered.
    % Axis/title text objects are excluded so only the atom labels are
    % affected.
        if isnan(fsz) || isempty(fig) || ~isgraphics(fig)
            return
        end
        ax = findall(fig, 'Type', 'Axes');
        if isempty(ax)
            return
        end
        ax = ax(1);

        txObjs = findall(ax, 'Type', 'Text');
        if isempty(txObjs)
            return
        end

        % exclude the axes title/xlabel/ylabel/zlabel, which are also
        % 'Text' objects but are not atom labels
        nonLabel = gobjects(0);
        for propName = ["Title","XLabel","YLabel","ZLabel"]
            h = ax.(propName);
            if isgraphics(h)
                nonLabel(end+1) = h; %#ok<AGROW>
            end
        end
        isAtomLabel = ~ismember(txObjs, nonLabel);
        atomLabelObjs = txObjs(isAtomLabel);

        fprintf('G09_modeViewer: %d Text object(s) found, %d treated as atom labels (FontSize -> %g).\n', ...
            numel(txObjs), numel(atomLabelObjs), fsz);

        if ~isempty(atomLabelObjs)
            set(atomLabelObjs, 'FontSize', fsz);
        end
    end

% -------------------------------------------------------------------------
    function onOptionChanged()
    % Any Scale/AtomScale/BondTol/ShowLabels/FlipSign control changed:
    % close the current figure and redraw the same mode with the updated
    % settings (in place, no accumulation of windows).
        drawMode(modeDD.Value, true);
    end

% -------------------------------------------------------------------------
    function onColorChanged(src)
    % ArrowColor combo box changed. 'Custom...' opens a colour picker.
    % Like other option changes, this replaces the current figure.
        if strcmp(src.Value, 'Custom...')
            picked = uisetcolor(customColor, 'Pick arrow colour');
            customColor = picked;   % uisetcolor returns the old colour if cancelled
        end
        drawMode(modeDD.Value, true);
    end

% -------------------------------------------------------------------------
    function applyTitle()
    % Sets the text box contents as the title of the target mode figure
    % (both the axes title text and the figure window name).
        target = resolveTargetFigure();
        if isempty(target) || ~isgraphics(target)
            uialert(selFig, 'No mode figure is currently open.', 'Nothing to title');
            return
        end
        newTitle = titleEdit.Value;
        if isempty(newTitle)
            return
        end
        ax = findall(target, 'Type', 'Axes');
        if ~isempty(ax)
            title(ax(1), newTitle, 'Interpreter', 'none');
        end
        target.Name = newTitle;
    end

% -------------------------------------------------------------------------
    function saveCurrentFigure(fmt)
    % Exports the target mode figure to JPEG/EPS/PDF via exportgraphics.
    % Vector output (true, editable vector paths) is used for PDF and EPS.
        target = resolveTargetFigure();
        if isempty(target) || ~isgraphics(target)
            uialert(selFig, 'No mode figure is currently open.', 'Nothing to save');
            return
        end

        switch fmt
            case 'pdf',   ext = '.pdf';  filterSpec = {'*.pdf',  'PDF file (vector)'};
            case 'epsc',  ext = '.eps';  filterSpec = {'*.eps',  'EPS file (vector, color)'};
            case 'jpeg',  ext = '.jpg';  filterSpec = {'*.jpg',  'JPEG file (raster)'};
            otherwise
                ext = '.pdf'; filterSpec = {'*.pdf', 'PDF file (vector)'};
        end

        defaultName = [matlab.lang.makeValidName(target.Name), ext];
        [f, p] = uiputfile(filterSpec, 'Save mode figure as', defaultName);
        if isequal(f, 0)
            return   % user cancelled
        end
        outFile = fullfile(p, f);

        try
            if strcmp(fmt, 'jpeg')
                exportgraphics(target, outFile, 'Resolution', 300);
            else
                exportgraphics(target, outFile, 'ContentType', 'vector');
            end
            fprintf('G09_modeViewer: figure saved to %s\n', outFile);
        catch ME
            uialert(selFig, ME.message, 'Export error');
        end
    end

% -------------------------------------------------------------------------
    function animateCurrentMode()
    % Exports an MP4 animation of the currently selected mode via
    % G09_animate_mode, using the Scale/AtomScale/BondTol/ShowLabels/
    % FlipSign values currently set in the option controls (ArrowColor is
    % not applicable: the animation shows the oscillating structure only,
    % no displacement arrows). The animation starts from whatever camera
    % orientation the currently displayed mode figure is in (so a manual
    % rotation carries over), falling back to MATLAB's default 3D view if
    % no mode figure is open or its view cannot be read.
        k = modeDD.Value;
        defaultName = sprintf('%s_mode%d.mp4', matlab.lang.makeValidName(filename), k);
        [f, p] = uiputfile({'*.mp4', 'MP4 video'}, 'Save mode animation as', defaultName);
        if isequal(f, 0)
            return   % user cancelled
        end
        outFile = fullfile(p, f);

        viewAngle = [];
        target = resolveTargetFigure();
        if ~isempty(target) && isgraphics(target)
            axTarget = findall(target, 'Type', 'Axes');
            if ~isempty(axTarget)
                [az, el] = view(axTarget(1));
                viewAngle = [az, el];
            end
        end

        try
            G09_animate_mode(mol, nm, k, ...
                'Filename',   outFile, ...
                'Scale',      str2double(scaleDD.Value), ...
                'AtomScale',  str2double(atomScaleDD.Value), ...
                'BondTol',    str2double(bondTolDD.Value), ...
                'ShowLabels', strcmp(showLabelsDD.Value, 'true'), ...
                'FlipSign',   strcmp(flipSignDD.Value, 'true'), ...
                'View',       viewAngle);
            fprintf('G09_modeViewer: animation saved to %s\n', outFile);
        catch ME
            uialert(selFig, ME.message, 'G09_animate_mode error');
        end
    end

% -------------------------------------------------------------------------
    function fig = resolveTargetFigure()
    % Prefers whichever mode figure the user last clicked on (groot's
    % CurrentFigure); falls back to the most recently drawn one.
        openModeFigs = openModeFigs(isgraphics(openModeFigs));  % drop stale handles
        cf = groot().CurrentFigure;
        if ~isempty(cf) && isgraphics(cf) && any(cf == openModeFigs)
            fig = cf;
        elseif ~isempty(currentModeFig) && isgraphics(currentModeFig)
            fig = currentModeFig;
        elseif ~isempty(openModeFigs)
            fig = openModeFigs(end);
        else
            fig = gobjects(0);
        end
    end

% -------------------------------------------------------------------------
    function closeViewer()
    % Closes every mode figure opened by this viewer, then the selector.
        openModeFigs = openModeFigs(isgraphics(openModeFigs));
        for i = 1:numel(openModeFigs)
            close(openModeFigs(i));
        end
        delete(selFig);
    end
% =========================================================================

end % function G09_modeViewer


% ------------------------------------------------------------------------
function s = tf2str(b)
%TF2STR  'true'/'false' string for a logical value.
if b, s = 'true'; else, s = 'false'; end
end