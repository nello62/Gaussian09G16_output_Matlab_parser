function G_raman_browser(filename, varargin)
% G_RAMAN_BROWSER  Click-to-select stick spectrum browser for vibrational
%                   modes: plots intensity vs. frequency as a stem
%                   ("stick") plot, click a stick to select the nearest
%                   mode, then press a button to render its 3D
%                   displacement structure.
%
%   G_raman_browser(filename)
%   G_raman_browser(filename, 'Name', Value, ...)
%
%   This is a standalone G_Utility companion to the core toolbox's
%   G16_MODEVIEWER/G09_MODEVIEWER: those tools select a mode from a text
%   drop-down list; this one selects it by clicking directly on its
%   stick in the spectrum, which is often the more natural way to pick
%   out "that strong peak at ~1580 cm-1" without first checking the
%   value in a list.
%
%   Works with BOTH Gaussian 09 and Gaussian 16 output files: the
%   Gaussian major version is auto-detected from the "Gaussian NN,
%   Revision X.YY," citation line near the top of the file (via
%   G16_GAUSSIAN_VERSION, which despite its name reads this for any NN),
%   and G09_STRUCTURE/G09_NMODES/G09_DRAW_MODE or
%   G16_STRUCTURE/G16_NMODES/G16_DRAW_MODE are used accordingly -- their
%   output structs/signatures are field-for-field identical, so the rest
%   of this function is version-agnostic. If the version cannot be
%   determined, G16 is assumed (this function's original, and still
%   most common, use case). No core-toolbox modification either way.
%
%   Input:
%       filename - a Gaussian output file (.log/.out), from either
%                  Gaussian 09 or Gaussian 16.
%
%   Optional parameters (Name-Value):
%       'Property'  - 'Raman' (default) | 'IR' -- which stick intensity
%                     to plot. 'Raman' requires the file to actually
%                     contain Raman activities (a Freq=Raman job); a
%                     clear error suggests 'IR' instead otherwise.
%       'Scale'     - forwarded to G16_draw_mode (default: 1.5)
%       'AtomScale' - forwarded to G16_draw_mode (default: 0.35)
%       'ShowLabels'- forwarded to G16_draw_mode (default: false)
%
%   Click anywhere in the spectrum axes to select the nearest mode by
%   frequency (highlighted in red); the info line above the plot updates
%   immediately. Press "Draw mode" to open (or replace) a 3D structure
%   figure showing that mode's displacement arrows.
%
%   Example:
%       G_raman_browser('4-NTP.out');                          % G16 file
%       G_raman_browser('indaco.log', 'Property', 'IR');       % G09 file
%       G_raman_browser('4-NTP.out', 'Property', 'IR', 'Scale', 2);
%
%   See also G16_MODEVIEWER, G16_DRAW_MODE, G16_NMODES, G16_STRUCTURE,
%            G09_MODEVIEWER, G09_DRAW_MODE, G09_NMODES, G09_STRUCTURE,
%            G16_GAUSSIAN_VERSION.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

if ~(ischar(filename) || (isstring(filename) && isscalar(filename)))
    error('G_raman_browser:invalidInput', ...
        'G_raman_browser expects a Gaussian output FILENAME (char/string). Usage: G_raman_browser(''molecule.out'')');
end

p = inputParser;
addRequired(p,  'filename', @(s) true);
addParameter(p, 'Property',   'Raman', @(s) any(strcmpi(s, {'raman','ir'})));
addParameter(p, 'Scale',      1.5,     @isnumeric);
addParameter(p, 'AtomScale',  0.35,    @isnumeric);
addParameter(p, 'ShowLabels', false,   @islogical);
parse(p, filename, varargin{:});

prop        = lower(p.Results.Property);
drawScale   = p.Results.Scale;
drawAtomScl = p.Results.AtomScale;
drawLabels  = p.Results.ShowLabels;

% -------------------------------------------------------------------------
% Detect Gaussian 09 vs Gaussian 16 and dispatch to the matching toolbox
% functions. G09_structure/G16_structure, G09_nmodes/G16_nmodes, and
% G09_draw_mode/G16_draw_mode are field-for-field/argument-for-argument
% identical in their outputs and signatures, so a plain function-handle
% swap is all that's needed -- everything below this point is written
% once, against whichever pair was selected here.
% -------------------------------------------------------------------------
gv = G16_gaussian_version(filename);
if isequal(gv.major, 9)
    structureFcn = @G09_structure;
    nmodesFcn    = @G09_nmodes;
    drawModeFcn  = @G09_draw_mode;
    verLabel     = 'Gaussian 09';
else
    structureFcn = @G16_structure;
    nmodesFcn    = @G16_nmodes;
    drawModeFcn  = @G16_draw_mode;
    verLabel     = 'Gaussian 16';
    if isempty(gv.major)
        fprintf('G_raman_browser: Gaussian version not detected in %s, assuming Gaussian 16.\n', filename);
    end
end

fprintf('G_raman_browser: reading structure and normal modes from %s (%s) ...\n', filename, verLabel);
mol = structureFcn(filename);
nm  = nmodesFcn(filename);
fprintf('  %d atoms, %d vibrational modes.\n', nm.Natoms, nm.Nmodes);

if nm.Nmodes < 1
    error('G_raman_browser:noModes', 'No vibrational modes found in %s.', filename);
end
if strcmp(prop, 'raman') && ~nm.has_Raman
    error('G_raman_browser:noRaman', ...
        ['%s has no Raman activities (not a Freq=Raman job) -- pass ''Property'',''IR'' ' ...
         'to browse the IR spectrum instead.'], filename);
end

if strcmp(prop, 'raman')
    vals = nm.Raman;
    ylab = 'Raman activity (Å^4/AMU)';
else
    vals = nm.IR;
    ylab = 'IR intensity (KM/Mole)';
end
freqs = nm.freq;

selectedIdx = [];       % currently selected mode index (into nm.*), [] = none yet
modeFig     = gobjects(0);   % the single mode-structure figure this browser drives

% -------------------------------------------------------------------------
% Build the browser window
% -------------------------------------------------------------------------
% A fixed 'Position' (e.g. [100 100 800 560]) is given in MATLAB's virtual
% desktop coordinates, which can land off-screen on multi-monitor setups
% where the primary monitor's origin is not at [0 0] (a real issue found
% and fixed the same way in G16_MODEVIEWER): the window is technically
% created but never visibly appears. To avoid that, a short-lived
% invisible ordinary FIGURE is used purely to ask MATLAB which monitor a
% new window would normally appear on right now, and the uifigure is then
% placed on THAT monitor in a single atomic call (no Position read/write
% on the uifigure itself, which is unsafe -- see G16_MODEVIEWER for the
% full explanation).
panelW = 800;
panelH = 560;
tmpFig = figure('Visible', 'off');

drawnow;
tmpPos = tmpFig.Position;
delete(tmpFig);

mp = get(groot, 'MonitorPositions');   % [Nmonitors x 4], each row [left bottom width height]
monIdx = find(tmpPos(1) >= mp(:,1) & tmpPos(1) <= mp(:,1) + mp(:,3) & ...
              tmpPos(2) >= mp(:,2) & tmpPos(2) <= mp(:,2) + mp(:,4), 1);
if isempty(monIdx), monIdx = 1; end
scr = mp(monIdx, :);                   % [left bottom width height] px, the active monitor

winX = scr(1) + 20;
winY = max(scr(2) + 40, scr(2) + scr(4) - panelH - 80);

fig = uifigure('Name', sprintf('G_raman_browser - %s', filename), ...
               'Position', [winX winY panelW panelH]);
infoLabel = uilabel(fig, 'Position', [20 525 760 24], ...
    'Text', sprintf('Click a stick to select a mode (%d modes, %s).', nm.Nmodes, prop), ...
    'FontWeight', 'bold', 'FontSize', 12);

ax = uiaxes(fig, 'Position', [20 90 760 420]);
hStem = stem(ax, freqs, vals, 'filled', 'Marker', 'o', 'MarkerSize', 4, ...
    'Color', [0.20 0.40 0.80], 'LineWidth', 1.1);
hStem.ButtonDownFcn = @(src, evt) onClick(evt);
ax.ButtonDownFcn = @(src, evt) onClick(evt);
hold(ax, 'on');
hHighlight = plot(ax, NaN, NaN, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.85 0.10 0.10], ...
    'MarkerEdgeColor', [0.85 0.10 0.10]);
hHighlight.ButtonDownFcn = @(src, evt) onClick(evt);
hold(ax, 'off');
xlabel(ax, 'Frequency (cm^{-1})');
ylabel(ax, ylab, 'Interpreter', 'tex');
ax.Interactions = [];   % disable pan/zoom gestures eating plain clicks
grid(ax, 'on');
box(ax, 'on');

uibutton(fig, 'push', 'Position', [20 40 200 32], ...
    'Text', 'Draw mode', 'FontWeight', 'bold', ...
    'ButtonPushedFcn', @(src, evt) drawSelectedMode());

uilabel(fig, 'Position', [240 40 540 32], ...
    'Text', sprintf('File: %s', filename), 'FontColor', [0.45 0.45 0.45]);
drawnow
    function onClick(evt)
        xClick = evt.IntersectionPoint(1);
        [~, k] = min(abs(freqs - xClick));
        selectedIdx = k;
        set(hHighlight, 'XData', freqs(k), 'YData', vals(k));
        symLabel = '';
        if ~isempty(nm.symmetry) && numel(nm.symmetry) >= k && ~isempty(nm.symmetry{k})
            symLabel = sprintf('  (%s)', nm.symmetry{k});
        end
        if strcmp(prop, 'raman') && nm.has_RamanFr
            % Both the static (0 cm-1) and the dynamic/CPHF=RdFreq
            % activities are available (G09/G16_nmodes' .RamanFr): show
            % all of them rather than just the static value that .Raman
            % (and a plain 'raman = %.2f') would otherwise silently be.
            frText = '';
            for c = 2:size(nm.RamanFr, 2)
                if nm.IncidentLight(c) > 0
                    wl_nm = 1e7 / nm.IncidentLight(c);
                    frText = [frText, sprintf(', %.0fnm = %.2f', wl_nm, nm.RamanFr(k, c))]; %#ok<AGROW>
                end
            end
            infoLabel.Text = sprintf('Selected: mode %d, %.1f cm^{-1}%s, Raman static = %.2f%s', ...
                k, freqs(k), symLabel, nm.RamanFr(k, 1), frText);
        else
            infoLabel.Text = sprintf('Selected: mode %d, %.1f cm^{-1}%s, %s = %.2f', ...
                k, freqs(k), symLabel, prop, vals(k));
        end
    end

    function drawSelectedMode()
        if isempty(selectedIdx)
            uialert(fig, 'Click a stick in the spectrum first to select a mode.', 'No mode selected');
            return
        end
        try
            drawModeFcn(mol, nm, selectedIdx, 'Scale', drawScale, ...
                'AtomScale', drawAtomScl, 'ShowLabels', drawLabels);
        catch ME
            uialert(fig, ME.message, 'draw_mode error');
            return
        end
        modeFig = gcf;
        modeFig.Name = sprintf('%s - Mode %d (%.1f cm^{-1})', filename, selectedIdx, freqs(selectedIdx));
    end

end % G_raman_browser
