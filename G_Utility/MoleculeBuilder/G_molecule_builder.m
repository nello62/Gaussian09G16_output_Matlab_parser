function G_molecule_builder()
% G_MOLECULE_BUILDER  Interactive molecule builder for the G09/G16 Toolbox:
%                     place atoms, form bonds, and export a Gaussian
%                     input (.gjf) file.
%
%   G_molecule_builder()
%
%   v1 scope (deliberately): no fragment/ring template library, no
%   automatic hybridization-aware geometry (tetrahedral/trigonal spread
%   for a 3rd+ substituent on the same atom) -- candidates for a later
%   version, once this one is validated in use. Every placement/edit is
%   one of:
%     - "Add first atom"            -> placed at the origin
%     - "Add atom bonded to selected" -> placed at
%       covalent_radius(selected)+covalent_radius(new) from the selected
%       atom, direction from NEWSUBSTITUENTDIRECTION (a generic,
%       non-degenerate spread -- 109.5 degrees from the first existing
%       bond, then away from the growing set after that; not a
%       hybridization solver -- refine afterward by dragging it, or via
%       the Atoms table's X/Y/Z columns, or "Optimize structure")
%     - "Bond selected pair"        -> select exactly two rows in the
%       Atoms table, then toggle a bond between them
%     - "Align A1->A2 to axis"      -> select exactly two rows in the
%       Atoms table (A1 = lower row index, A2 = higher one, since a
%       UITABLE's multi-selection is reported sorted, not in click order)
%       and rigidly rotate+translate the whole molecule so both A1 and A2
%       end up lying exactly ON the chosen axis (X/Y/Z) -- i.e. their two
%       off-axis coordinates both become exactly 0 (e.g. choosing X makes
%       A1_y=A1_z=A2_y=A2_z=0). Useful to orient a molecule along a
%       bond/axis of interest (e.g. before reading off a dipole component
%       by eye, or just for a cleaner view) -- purely a rigid motion,
%       bond lengths/angles are unaffected.
%     - "Align A1->A2 to vector"    -> same as the above, but the target
%       direction is an arbitrary vector [x,y,z] typed into the "Vec"
%       fields instead of a coordinate axis -- e.g. a molecule's own
%       dipole moment vector (as reported by a Gaussian job on it), to
%       view the molecule relative to that direction. A1 ends up exactly
%       on the line through the origin along that vector (its component
%       perpendicular to the vector is zeroed), A2 (and any other atom
%       already collinear with A1/A2) lands on it too, same "for free"
%       reasoning as the axis case -- see RIGIDALIGNPAIRTOVECTOR, which
%       both this and "Align A1->A2 to axis" call (a coordinate axis is
%       just one particular choice of vector).
%     - "Align A1-A2-A3 to plane"   -> select exactly three rows in the
%       Atoms table (A1/A2/A3 = ascending row index, same convention) and
%       rigidly rotate+translate the molecule so all three end up lying
%       exactly in the chosen coordinate plane (XY/XZ/YZ) -- e.g.
%       choosing XY makes A1_z=A2_z=A3_z=0. Errors out (via UIALERT) if
%       the three are collinear/coincident (no unique plane).
%       "Show X/Y/Z axes (dashed)" (auto-checked the first time either
%       alignment is done, but toggleable any time) overlays dashed
%       +/-X/Y/Z reference lines through the origin, colored red/green/blue
%       and labelled, so the alignment is visually obvious in the 3D view.
%     - Atoms/Bonds tables          -> direct numeric/symbol/bond-order
%       editing, and row selection (selects the corresponding atom)
%     - Click an atom directly in the 3D view -> selects it (same as
%       selecting its row). With "Move atoms" checked, click-and-drag an
%       atom instead moves it, confined to the screen-parallel plane
%       through its position at the start of the drag (a 2D mouse
%       position alone cannot fix a depth without some constraint; this
%       is the standard choice for this kind of free-rotation 3D view --
%       rotate the view first to approach a desired depth/angle, then
%       drag). "Move atoms" and the view's own interactive rotation
%       (ROTATE3D, set up by G16_DRAW_MOLECULE itself every redraw) are
%       mutually exclusive by design, toggled by the same checkbox: a
%       MATLAB figure interaction mode (rotate3d/zoom/pan) is documented
%       to take over mouse events on the axes while enabled, which would
%       otherwise make an object's own ButtonDownFcn unreliable -- so
%       ROTATE3D is explicitly disabled (its Enable set back to 'off'
%       right after G16_draw_molecule turns it back on) for as long as
%       "Move atoms" stays checked, and re-enabled when unchecked. While
%       actively dragging, the bottom status strip shows the atom's live
%       x/y/z as the mouse moves (updated every ONMOUSEMOVE call, i.e.
%       every frame of the drag).
%
%   The 3D view keeps a visible outer frame (axes box) at all times --
%   G16_DRAW_MOLECULE itself calls AXIS(ax,'off') on every redraw (hiding
%   box/ticks/labels), so REDRAW turns AX.Visible/AX.Box back on right
%   after, every time, the same "re-assert after every call" pattern
%   already used for ROTATE3D's Enable property above.
%
%   Reuses, rather than reimplements:
%     - G16_DRAW_MOLECULE for the entire live 3D CPK preview (it already
%       accepts an arbitrary 'mol' struct plus an explicit 'BondList'
%       [Nbonds x 3] Atom1/Atom2/BondOrder override -- the same
%       convention G16_NBO_BONDS/G16_CHARGES already use)
%     - G_ELEMENT_TABLE for covalent radii (bond-length defaults) --
%       consolidates G16_GET_BOND_LENGTH's and G16_DRAW_MOLECULE's two
%       independent embedded element tables into one shared source; those
%       two existing, already-shipped functions are deliberately NOT
%       refactored to call it in this change (avoids touching working
%       code for a DRY-only benefit; a natural, low-risk follow-up)
%     - G_WRITE_GJF for the actual .gjf export (factored out of
%       G16_RESTART's own inline writing block, same rationale)
%     - G_READ_INPUT + G16_GET_BOND_LENGTH for the optional "Load input..."
%       import (.gjf/.com/.in -- same free-form Gaussian input text format
%       regardless of extension, so no extension-specific parsing is
%       needed). Starting bonds are geometrically auto-detected from the
%       loaded geometry, same heuristic G16_GET_BOND_LENGTH already uses
%       elsewhere -- editable afterward via the Bonds table, same as any
%       other bond here)
%     - "Optimize structure" dispatches, per the "Method" dropdown next to
%       it, to one of:
%         - G_ROUGH_OPTIMIZE_GEOMETRY ("Quick cleanup (built-in)"): a
%           quick, approximate geometry clean-up (bond lengths toward
%           covalent-radii sums, angles spread toward a generic
%           tetrahedral-like default, non-bonded clashes pushed apart),
%           via LSQNONLIN -- not a real force field, no external
%           dependency.
%         - G_OPENBABEL_OPTIMIZE_GEOMETRY ("Open Babel: MMFF94/MMFF94s/
%           UFF"): a real, validated force-field minimization via Open
%           Babel (https://openbabel.org) -- the one feature in this app
%           with a non-MATLAB dependency; see that function's own error
%           message (shown in a dialog if it is not found) for how to
%           install it.
%       Neither is a substitute for a real Gaussian 'opt' job -- both just
%       de-crumple a hand-placed sketch first.
%
%   See also G16_DRAW_MOLECULE, G_ELEMENT_TABLE, G_WRITE_GJF,
%            G_READ_INPUT, G16_GET_BOND_LENGTH, G_ROUGH_OPTIMIZE_GEOMETRY,
%            G_OPENBABEL_OPTIMIZE_GEOMETRY, G_GAUSSIAN_VIEWER.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

% -------------------------------------------------------------------------
% Session state (nested-function closures share these, same pattern as
% G_GAUSSIAN_VIEWER / G16_MODEVIEWER -- no classdef app for this tool).
% -------------------------------------------------------------------------
symbols  = {};          % {Natoms x 1 cell}
xyz      = zeros(0, 3); % [Natoms x 3]
bondList = zeros(0, 3); % [Nbonds x 3] Atom1, Atom2, Order

selAtomIdx          = 0;   % single selected atom row (0 = none)
lastAtomsSelection  = [];  % all currently-selected Atoms-table rows (for Bond pair)
lastBondsSelection  = [];  % all currently-selected Bonds-table rows (for Remove bond)

isDragging       = false;  % true while a 3D-view atom drag is in progress
dragIdx          = 0;      % atom index being dragged
dragPlanePoint   = [0 0 0];% drag constraint plane: a point on it (atom's position when the drag started)
dragPlaneNormal  = [0 0 1];% drag constraint plane: its normal (view direction when the drag started)

historyStack = {struct('symbols', {{}}, 'xyz', zeros(0,3), 'bondList', zeros(0,3))};
historyIdx   = 1;

[cpkColorsTbl, covRadiiTbl] = G_element_table(); %#ok<ASGLU> -- cpkColorsTbl unused here, kept for symmetry/future use

% -------------------------------------------------------------------------
% Find the active monitor, then build the window in a single atomic call
% (same "short-lived invisible ordinary figure" trick as G_GAUSSIAN_VIEWER/
% G16_MODEVIEWER/G_RAMAN_BROWSER).
% -------------------------------------------------------------------------
winW = 1200;
winH = 900;
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

fig = uifigure('Name', 'G_Molecule_Builder', 'Position', [winX winY winW winH]);
fig.CloseRequestFcn = @(src, evt) delete(fig);

% -------------------------------------------------------------------------
% Window layout: bottom status strip, left sidebar, central 3D view.
% -------------------------------------------------------------------------
sidebarW = 380;

statusLabel = uilabel(fig, 'Position', [10 8 winW-20 28], ...
    'Text', 'New molecule. Pick an element and click "Add first atom" to begin.', ...
    'FontColor', [0.35 0.35 0.35]);

sidebar = uipanel(fig, 'Position', [0 40 sidebarW winH-40], 'BorderType', 'line');
ax = uiaxes(fig, 'Position', [sidebarW+10 45 winW-sidebarW-20 winH-50]);
ax.Toolbar.Visible = 'on';
% Interactive rotation is enabled by G16_draw_molecule itself (it sets up
% rotate3d, with its own camlight-relinking callback, every time it
% renders into AX) -- no separate setup needed here. The figure-level
% motion/up handlers for the atom-drag itself are (re-)wired inside
% REDRAW(), not just once here -- see its own comment for why.

% -------------------------------------------------------------------------
% Sidebar: fixed top controls (element picker, add/remove), tab group below
% -------------------------------------------------------------------------
uilabel(sidebar, 'Position', [10 830 sidebarW-20 24], ...
    'Text', 'Molecule Builder', 'FontWeight', 'bold', 'FontSize', 14);

uibutton(sidebar, 'push', 'Position', [10 798 175 28], ...
    'Text', 'New molecule', 'ButtonPushedFcn', @(s,e) onNewMolecule());
uibutton(sidebar, 'push', 'Position', [195 798 175 28], ...
    'Text', 'Load input...', 'ButtonPushedFcn', @(s,e) onLoadGjf());

uilabel(sidebar, 'Position', [10 768 100 18], 'Text', 'Element:', 'FontWeight', 'bold');
elementDD = uidropdown(sidebar, 'Position', [10 744 sidebarW-30 24], ...
    'Items', {'H','C','N','O','F','P','S','Cl','Br','I','Si','B','Other...'}, 'Value', 'C');
uilabel(sidebar, 'Position', [10 716 115 18], 'Text', 'Custom symbol:');
customSymField = uieditfield(sidebar, 'text', 'Position', [125 714 sidebarW-155 22], 'Value', '');

addFirstBtn = uibutton(sidebar, 'push', 'Position', [10 682 175 28], ...
    'Text', 'Add first atom', 'ButtonPushedFcn', @(s,e) onAddFirstAtom());
addBondedBtn = uibutton(sidebar, 'push', 'Position', [195 682 175 28], ...
    'Text', 'Add bonded atom', 'Enable', 'off', 'ButtonPushedFcn', @(s,e) onAddBondedAtom());
removeAtomBtn = uibutton(sidebar, 'push', 'Position', [10 650 sidebarW-30 28], ...
    'Text', 'Remove selected atom', 'Enable', 'off', 'ButtonPushedFcn', @(s,e) onRemoveAtom());

moveModeCheck = uicheckbox(sidebar, 'Position', [10 620 sidebarW-30 22], ...
    'Text', 'Move atoms (click + drag in 3D view; disables view rotation while on)', ...
    'Value', false, 'ValueChangedFcn', @(s,e) redraw());

tg = uitabgroup(sidebar, 'Position', [5 10 sidebarW-10 605]);
tabBuild = uitab(tg, 'Title', 'Atoms & Bonds');
tabJob   = uitab(tg, 'Title', 'Gaussian Job');

% ---- "Atoms & Bonds" tab -------------------------------------------------
uilabel(tabBuild, 'Position', [5 600 150 18], 'Text', 'Atoms:', 'FontWeight', 'bold');
atomsTable = uitable(tabBuild, 'Position', [5 360 sidebarW-40 235], ...
    'ColumnName', {'Idx','Symbol','X','Y','Z'}, ...
    'ColumnEditable', [false true true true true], ...
    'CellEditCallback', @(s,e) onAtomsTableEdit(e), ...
    'CellSelectionCallback', @(s,e) onAtomsTableSelection(e));

bondPairBtn = uibutton(tabBuild, 'push', 'Position', [5 325 160 28], ...
    'Text', 'Bond selected pair', 'Enable', 'off', 'ButtonPushedFcn', @(s,e) onBondPair());
removeBondBtn = uibutton(tabBuild, 'push', 'Position', [170 325 160 28], ...
    'Text', 'Remove selected bond', 'Enable', 'off', 'ButtonPushedFcn', @(s,e) onRemoveBond());

uilabel(tabBuild, 'Position', [5 296 40 18], 'Text', 'Axis:');
alignAxisDD = uidropdown(tabBuild, 'Position', [48 293 55 24], ...
    'Items', {'X','Y','Z'}, 'Value', 'Z');
alignBtn = uibutton(tabBuild, 'push', 'Position', [108 293 sidebarW-143 24], ...
    'Text', 'Align A1->A2 to axis', 'Enable', 'off', 'ButtonPushedFcn', @(s,e) onAlignToAxis());

uilabel(tabBuild, 'Position', [5 264 40 18], 'Text', 'Plane:');
alignPlaneDD = uidropdown(tabBuild, 'Position', [48 261 55 24], ...
    'Items', {'XY','XZ','YZ'}, 'Value', 'XY');
alignPlaneBtn = uibutton(tabBuild, 'push', 'Position', [108 261 sidebarW-143 24], ...
    'Text', 'Align A1-A2-A3 to plane', 'Enable', 'off', 'ButtonPushedFcn', @(s,e) onAlignToPlane());

uilabel(tabBuild, 'Position', [5 234 30 18], 'Text', 'Vec:');
vecXField = uieditfield(tabBuild, 'numeric', 'Position', [38 231 44 24], 'Value', 0);
vecYField = uieditfield(tabBuild, 'numeric', 'Position', [85 231 44 24], 'Value', 0);
vecZField = uieditfield(tabBuild, 'numeric', 'Position', [132 231 44 24], 'Value', 1);
alignVecBtn = uibutton(tabBuild, 'push', 'Position', [179 231 sidebarW-174 24], ...
    'Text', 'Align A1->A2 to vector', 'Enable', 'off', 'ButtonPushedFcn', @(s,e) onAlignToVector());

showAxesCheck = uicheckbox(tabBuild, 'Position', [5 201 sidebarW-40 22], ...
    'Text', 'Show X/Y/Z axes (dashed)', 'Value', false, 'ValueChangedFcn', @(s,e) redraw());

uilabel(tabBuild, 'Position', [5 175 150 18], 'Text', 'Bonds:', 'FontWeight', 'bold');
bondsTable = uitable(tabBuild, 'Position', [5 120 sidebarW-40 50], ...
    'ColumnName', {'Atom1','Atom2','Order'}, ...
    'ColumnEditable', [false false true], ...
    'CellEditCallback', @(s,e) onBondsTableEdit(e), ...
    'CellSelectionCallback', @(s,e) onBondsTableSelection(e));

optMethodDD = uidropdown(tabBuild, 'Position', [5 90 sidebarW-40 22], ...
    'Items', {'Quick cleanup (built-in)', 'Open Babel: MMFF94', 'Open Babel: MMFF94s', 'Open Babel: UFF'}, ...
    'Value', 'Quick cleanup (built-in)');
uibutton(tabBuild, 'push', 'Position', [5 57 sidebarW-40 28], ...
    'Text', 'Optimize structure', 'ButtonPushedFcn', @(s,e) onOptimizeStructure());

undoBtn = uibutton(tabBuild, 'push', 'Position', [5 19 160 30], ...
    'Text', 'Undo', 'Enable', 'off', 'ButtonPushedFcn', @(s,e) onUndo());
redoBtn = uibutton(tabBuild, 'push', 'Position', [170 19 160 30], ...
    'Text', 'Redo', 'Enable', 'off', 'ButtonPushedFcn', @(s,e) onRedo());

% ---- "Gaussian Job" tab ---------------------------------------------------
uilabel(tabJob, 'Position', [5 555 150 18], 'Text', 'Route:');
routeField = uieditfield(tabJob, 'text', 'Position', [5 530 sidebarW-40 24], ...
    'Value', '# opt freq b3lyp/6-31g(d)');
uilabel(tabJob, 'Position', [5 500 150 18], 'Text', 'Title:');
titleField = uieditfield(tabJob, 'text', 'Position', [5 475 sidebarW-40 24], 'Value', 'Built with G_molecule_builder');

uilabel(tabJob, 'Position', [5 440 80 18], 'Text', 'Charge:');
chargeField = uieditfield(tabJob, 'numeric', 'Position', [90 438 100 24], ...
    'Value', 0, 'RoundFractionalValues', 'on');
uilabel(tabJob, 'Position', [200 440 60 18], 'Text', 'Mult.:');
multField = uieditfield(tabJob, 'numeric', 'Position', [255 438 90 24], ...
    'Value', 1, 'Limits', [1 Inf], 'RoundFractionalValues', 'on');

uilabel(tabJob, 'Position', [5 400 sidebarW-40 18], 'Text', 'Link0 (optional):', 'FontWeight', 'bold');
uilabel(tabJob, 'Position', [5 375 50 18], 'Text', '%chk:');
chkField = uieditfield(tabJob, 'text', 'Position', [60 373 sidebarW-95 22], 'Value', '');
uilabel(tabJob, 'Position', [5 345 50 18], 'Text', '%mem:');
memField = uieditfield(tabJob, 'text', 'Position', [60 343 sidebarW-95 22], 'Value', '');
uilabel(tabJob, 'Position', [5 315 100 18], 'Text', '%nprocshared:');
nprocField = uieditfield(tabJob, 'numeric', 'Position', [110 313 sidebarW-145 22], ...
    'Value', 0, 'Limits', [0 Inf], 'RoundFractionalValues', 'on');

uibutton(tabJob, 'push', 'Position', [5 20 sidebarW-40 36], ...
    'Text', 'Generate Gaussian input...', 'FontWeight', 'bold', ...
    'ButtonPushedFcn', @(s,e) onGenerateInput());

updateAtomsTable();
updateBondsTable();

% =========================================================================
%  Nested callback/helper functions (share the outer function's workspace)
% =========================================================================

    function clearAx()
    % Removes every plotted object from the shared AX before a new
    % render -- CLA(ax) alone does not fully clear a uiaxes (some CPK
    % spheres/bond lines survive it); copied verbatim from
    % G_GAUSSIAN_VIEWER's own CLEARAX, same rationale.
        kids = findall(ax);
        kids(kids == ax) = [];
        delete(kids);
    end

% -------------------------------------------------------------------------
    function [f, p] = pickOpenFile(filterSpec, dlgTitle)
    % Same macOS "minimise the uifigure before UIGETFILE" workaround as
    % G_GAUSSIAN_VIEWER's PICKOPENFILE.
        prevState = fig.WindowState;
        if strcmp(prevState, 'minimized')
            prevState = 'normal';
        end
        fig.WindowState = 'minimized';
        drawnow;
        try
            [f, p] = uigetfile(filterSpec, dlgTitle);
        catch ME
            fig.WindowState = prevState;
            drawnow;
            rethrow(ME);
        end
        fig.WindowState = prevState;
        drawnow;
    end

% -------------------------------------------------------------------------
    function [f, p] = pickSaveFile(filterSpec, dlgTitle, defaultName)
    % Same rationale as PICKOPENFILE, for UIPUTFILE.
        prevState = fig.WindowState;
        if strcmp(prevState, 'minimized')
            prevState = 'normal';
        end
        fig.WindowState = 'minimized';
        drawnow;
        try
            [f, p] = uiputfile(filterSpec, dlgTitle, defaultName);
        catch ME
            fig.WindowState = prevState;
            drawnow;
            rethrow(ME);
        end
        fig.WindowState = prevState;
        drawnow;
    end

% -------------------------------------------------------------------------
    function sym = normalizeSymbol(raw)
    % Same element-symbol normalization convention G_READ_INPUT uses when
    % parsing a geometry line: first letter uppercase, optional second
    % letter lowercase (e.g. 'cl' / 'CL' -> 'Cl'). Returns '' if RAW does
    % not start with 1-2 letters.
        raw = strtrim(raw);
        tok = regexp(raw, '^([A-Za-z]{1,2})', 'tokens', 'once');
        if isempty(tok)
            sym = '';
            return
        end
        s = tok{1};
        if numel(s) > 1
            sym = [upper(s(1)), lower(s(2))];
        else
            sym = upper(s);
        end
    end

% -------------------------------------------------------------------------
    function sym = resolveElement()
        v = elementDD.Value;
        if strcmp(v, 'Other...')
            sym = normalizeSymbol(customSymField.Value);
        else
            sym = v;
        end
    end

% -------------------------------------------------------------------------
    function r = lookupRadius(sym)
    % Same fallback default (1.50 Angstrom) as G16_GET_BOND_LENGTH's own
    % LOOKUP_RADIUS, for an element not in G_ELEMENT_TABLE's covalent-radii map.
        if isKey(covRadiiTbl, sym)
            r = covRadiiTbl(sym);
        else
            r = 1.50;
        end
    end

% -------------------------------------------------------------------------
    function redraw()
        clearAx();
        if isempty(symbols)
            return
        end
        % .Natoms is required too, even though only .symbols/.xyz are
        % validated: G16_draw_molecule's own geometric bond-auto-detection
        % branch (ISEMPTY(bond_list) -- which BONDLIST hits here whenever
        % 0 bonds are declared yet, see the comment below) reads
        % mol.Natoms directly, uncondionally -- confirmed the hard way
        % (omitting it threw "Unrecognized field name 'Natoms'" the
        % moment a first atom was added, since bondList is then empty).
        mol = struct('symbols', {symbols(:)}, 'xyz', xyz, 'Natoms', numel(symbols));
        try
            % BONDLIST is always passed explicitly (even when empty) so
            % G16_draw_molecule respects the user's own declared bonds
            % rather than guessing from geometry -- its own fallback
            % (geometric auto-detection) only ever kicks in here while a
            % molecule has zero declared bonds yet, which is harmless
            % (purely cosmetic in the live preview; never written to the
            % exported .gjf, which has no bond section at all).
            G16_draw_molecule(mol, 'Ax', ax, 'BondList', bondList, ...
                'ShowLabels', true, 'ShowLegend', true, 'Title', '');
        catch ME
            uialert(fig, ME.message, 'G16_draw_molecule error');
            return
        end
        % G16_draw_molecule calls AXIS(ax,'off') itself on every call
        % (hides the box/ticks/labels entirely) -- turned back on here,
        % every redraw, so the 3D view always keeps a visible outer frame
        % (same "re-assert after every call" pattern REDRAW already uses
        % for ROTATE3D's own Enable property, for the same reason: the
        % upstream function resets it unconditionally each time).
        ax.Visible = 'on';
        ax.Box = 'on';
        ax.XGrid = 'off'; ax.YGrid = 'off'; ax.ZGrid = 'off';
        if selAtomIdx >= 1 && selAtomIdx <= numel(symbols)
            hold(ax, 'on');
            plot3(ax, xyz(selAtomIdx,1), xyz(selAtomIdx,2), xyz(selAtomIdx,3), 'o', ...
                'MarkerSize', 32, 'MarkerEdgeColor', [1.00 0.55 0.00], 'LineWidth', 2.5, ...
                'HandleVisibility', 'off');
        end
        if showAxesCheck.Value
            hold(ax, 'on');
            ext = max([1, max(abs(xyz(:))) * 1.3]);
            plot3(ax, [-ext ext], [0 0], [0 0], '--', 'Color', [0.85 0.10 0.10], ...
                'LineWidth', 1, 'HandleVisibility', 'off');
            plot3(ax, [0 0], [-ext ext], [0 0], '--', 'Color', [0.10 0.60 0.10], ...
                'LineWidth', 1, 'HandleVisibility', 'off');
            plot3(ax, [0 0], [0 0], [-ext ext], '--', 'Color', [0.10 0.30 0.85], ...
                'LineWidth', 1, 'HandleVisibility', 'off');
            text(ax, ext, 0, 0, 'X', 'Color', [0.85 0.10 0.10], 'FontWeight', 'bold', 'HandleVisibility', 'off');
            text(ax, 0, ext, 0, 'Y', 'Color', [0.10 0.60 0.10], 'FontWeight', 'bold', 'HandleVisibility', 'off');
            text(ax, 0, 0, ext, 'Z', 'Color', [0.10 0.30 0.85], 'FontWeight', 'bold', 'HandleVisibility', 'off');
        end
        tagAtomSpheres();
        % G16_draw_molecule just (re-)enabled ROTATE3D on AX (it does so
        % unconditionally on every call) -- turn it back off here while
        % "Move atoms" is checked, so a click-drag on a sphere moves the
        % atom instead of orbiting the camera (see header comment).
        % ROTATE3D(ax) returns the SAME behavior object every time
        % (confirmed empirically), so this toggles the one G16_draw_molecule
        % itself just created, not a second, competing one.
        hRot = rotate3d(ax);
        if moveModeCheck.Value
            hRot.Enable = 'off';
            % ROTATE3D's own mode manager installs ITS OWN
            % WindowButtonMotionFcn/WindowButtonUpFcn on FIG to implement
            % rotation-by-drag (confirmed empirically: after it has been
            % enabled even once, FIG.WindowButtonUpFcn is no longer the
            % plain function_handle this app set, but a cell-style
            % legacy callback) -- disabling it again (above) does not
            % hand those figure-level callbacks back, so they are
            % reclaimed here, explicitly, on every redraw while "Move
            % atoms" is checked (which includes every redraw triggered
            % by ONMOUSEMOVE itself during an active drag).
            fig.WindowButtonMotionFcn = @(s,e) onMouseMove();
            fig.WindowButtonUpFcn     = @(s,e) onMouseUp();
        else
            hRot.Enable = 'on';
            % Leave FIG's own WindowButtonMotionFcn/WindowButtonUpFcn
            % alone here: ROTATE3D needs to manage them itself to
            % implement rotation-by-drag while enabled.
        end
    end

% -------------------------------------------------------------------------
    function tagAtomSpheres()
    % Attaches a ButtonDownFcn to each atom's CPK sphere (G16_draw_molecule
    % draws exactly one Type=='surface' object per atom, nothing else of
    % that type), so clicking one in the 3D view selects it (and, with
    % "Move atoms" checked, starts dragging it) -- the sphere objects
    % themselves carry no atom-index tag of their own (G16_draw_molecule
    % returns no handles), so each is matched back to a row of XYZ by its
    % own centroid instead, which is robust regardless of FINDALL's
    % traversal/creation order.
        spheres = findall(ax, 'Type', 'surface');
        for k = 1:numel(spheres)
            h = spheres(k);
            centroid = [mean(h.XData(:)), mean(h.YData(:)), mean(h.ZData(:))];
            dists = vecnorm(xyz - centroid, 2, 2);
            [~, idx] = min(dists);
            h.ButtonDownFcn = @(src, evt) onAtomButtonDown(idx, evt);
        end
    end

% -------------------------------------------------------------------------
    function onAtomButtonDown(idx, ~)
        selAtomIdx = idx;
        lastAtomsSelection = idx;
        try
            atomsTable.Selection = idx;  % best-effort table-highlight sync; harmless if unsupported
        catch
        end
        updateButtonStates();
        if moveModeCheck.Value
            isDragging = true;
            dragIdx = idx;
            dragPlanePoint = xyz(idx, :);
            viewDir = ax.CameraPosition - ax.CameraTarget;
            viewDirNorm = norm(viewDir);
            if viewDirNorm > 1e-9
                dragPlaneNormal = viewDir / viewDirNorm;
            else
                dragPlaneNormal = [0 0 1];
            end
        end
        redraw();
    end

% -------------------------------------------------------------------------
    function onMouseMove()
        if ~isDragging
            return
        end
        rayNear = ax.CurrentPoint(1, :);
        rayFar  = ax.CurrentPoint(2, :);
        rayDir  = rayFar - rayNear;
        denom = dot(rayDir, dragPlaneNormal);
        if abs(denom) < 1e-9
            return   % ray parallel to the drag plane -- degenerate, skip this frame
        end
        t = dot(dragPlanePoint - rayNear, dragPlaneNormal) / denom;
        xyz(dragIdx, :) = rayNear + t * rayDir;
        statusLabel.Text = sprintf('Moving atom %d (%s): x=%.4f, y=%.4f, z=%.4f', ...
            dragIdx, symbols{dragIdx}, xyz(dragIdx,1), xyz(dragIdx,2), xyz(dragIdx,3));
        updateAtomsTable();
        redraw();
    end

% -------------------------------------------------------------------------
    function onMouseUp()
        if ~isDragging
            return
        end
        isDragging = false;
        snapshotPush();
        statusLabel.Text = sprintf('Moved atom %d.', dragIdx);
    end

% -------------------------------------------------------------------------
    function updateAtomsTable()
        n = numel(symbols);
        data = cell(n, 5);
        for i = 1:n
            data{i,1} = int32(i);
            data{i,2} = symbols{i};
            data{i,3} = xyz(i,1);
            data{i,4} = xyz(i,2);
            data{i,5} = xyz(i,3);
        end
        atomsTable.Data = data;
    end

% -------------------------------------------------------------------------
    function updateBondsTable()
        n = size(bondList, 1);
        data = cell(n, 3);
        for i = 1:n
            data{i,1} = int32(bondList(i,1));
            data{i,2} = int32(bondList(i,2));
            data{i,3} = int32(bondList(i,3));
        end
        bondsTable.Data = data;
    end

% -------------------------------------------------------------------------
    function updateButtonStates()
        addFirstBtn.Enable     = matlab.lang.OnOffSwitchState(isempty(symbols));
        addBondedBtn.Enable    = matlab.lang.OnOffSwitchState(selAtomIdx >= 1);
        removeAtomBtn.Enable   = matlab.lang.OnOffSwitchState(selAtomIdx >= 1);
        bondPairBtn.Enable     = matlab.lang.OnOffSwitchState(numel(lastAtomsSelection) == 2);
        alignBtn.Enable        = matlab.lang.OnOffSwitchState(numel(lastAtomsSelection) == 2);
        alignVecBtn.Enable     = matlab.lang.OnOffSwitchState(numel(lastAtomsSelection) == 2);
        alignPlaneBtn.Enable   = matlab.lang.OnOffSwitchState(numel(lastAtomsSelection) == 3);
        removeBondBtn.Enable   = matlab.lang.OnOffSwitchState(~isempty(lastBondsSelection));
    end

% -------------------------------------------------------------------------
    function updateUndoRedoButtons()
        undoBtn.Enable = matlab.lang.OnOffSwitchState(historyIdx > 1);
        redoBtn.Enable = matlab.lang.OnOffSwitchState(historyIdx < numel(historyStack));
    end

% -------------------------------------------------------------------------
    function snapshotPush()
    % Call AFTER a mutation (and its redraw) completes, to record the new
    % current state. Any undone-then-abandoned redo branch is dropped,
    % same convention as a standard text-editor undo stack.
        historyStack = historyStack(1:historyIdx);
        historyStack{end+1} = struct('symbols', {symbols}, 'xyz', xyz, 'bondList', bondList);
        historyIdx = numel(historyStack);
        updateUndoRedoButtons();
    end

% -------------------------------------------------------------------------
    function restoreSnapshot(idx)
        snap     = historyStack{idx};
        symbols  = snap.symbols;
        xyz      = snap.xyz;
        bondList = snap.bondList;
        selAtomIdx = 0;
        lastAtomsSelection = [];
        lastBondsSelection = [];
        updateAtomsTable();
        updateBondsTable();
        redraw();
        updateButtonStates();
    end

% -------------------------------------------------------------------------
    function onUndo()
        if historyIdx <= 1
            return
        end
        historyIdx = historyIdx - 1;
        restoreSnapshot(historyIdx);
        updateUndoRedoButtons();
    end

% -------------------------------------------------------------------------
    function onRedo()
        if historyIdx >= numel(historyStack)
            return
        end
        historyIdx = historyIdx + 1;
        restoreSnapshot(historyIdx);
        updateUndoRedoButtons();
    end

% -------------------------------------------------------------------------
    function onNewMolecule()
        if isempty(symbols)
            return
        end
        sel = uiconfirm(fig, ...
            'Clear the current molecule and its undo history? This cannot be undone.', ...
            'New molecule', 'Options', {'Clear', 'Cancel'}, ...
            'DefaultOption', 2, 'CancelOption', 2, 'Icon', 'warning');
        if ~strcmp(sel, 'Clear')
            return
        end
        symbols = {}; xyz = zeros(0,3); bondList = zeros(0,3);
        selAtomIdx = 0; lastAtomsSelection = []; lastBondsSelection = [];
        historyStack = {struct('symbols', {{}}, 'xyz', zeros(0,3), 'bondList', zeros(0,3))};
        historyIdx = 1;
        updateAtomsTable();
        updateBondsTable();
        redraw();
        updateButtonStates();
        updateUndoRedoButtons();
        statusLabel.Text = 'New molecule. Pick an element and click "Add first atom" to begin.';
    end

% -------------------------------------------------------------------------
    function onAddFirstAtom()
        sym = resolveElement();
        if isempty(sym)
            uialert(fig, 'Enter a valid element symbol.', 'Invalid symbol');
            return
        end
        symbols = {sym}; xyz = [0 0 0]; bondList = zeros(0,3);
        selAtomIdx = 0;
        updateAtomsTable();
        updateBondsTable();
        redraw();
        updateButtonStates();
        snapshotPush();
        statusLabel.Text = sprintf('Added first atom: %s at the origin.', sym);
    end

% -------------------------------------------------------------------------
    function onAddBondedAtom()
        if selAtomIdx < 1
            return
        end
        sym = resolveElement();
        if isempty(sym)
            uialert(fig, 'Enter a valid element symbol.', 'Invalid symbol');
            return
        end
        bondLen = lookupRadius(symbols{selAtomIdx}) + lookupRadius(sym);

        existingPartners = [];
        if ~isempty(bondList)
            m1 = bondList(:,1) == selAtomIdx;
            m2 = bondList(:,2) == selAtomIdx;
            existingPartners = [bondList(m1,2); bondList(m2,1)];
        end
        existingDirs = zeros(0, 3);
        for k = 1:numel(existingPartners)
            d = xyz(existingPartners(k),:) - xyz(selAtomIdx,:);
            nv = norm(d);
            if nv > 1e-9
                existingDirs(end+1, :) = d / nv; %#ok<AGROW>
            end
        end
        dirVec = newSubstituentDirection(existingDirs);

        symbols{end+1, 1} = sym;
        xyz(end+1, :) = xyz(selAtomIdx,:) + dirVec * bondLen;
        newIdx = numel(symbols);
        bondList(end+1, :) = [selAtomIdx, newIdx, 1];

        updateAtomsTable();
        updateBondsTable();
        redraw();
        snapshotPush();
        statusLabel.Text = sprintf('Added %s bonded to atom %d (bond length %.3f A).', sym, selAtomIdx, bondLen);
    end

% -------------------------------------------------------------------------
    function dirVec = newSubstituentDirection(existingDirs)
    % Picks a unit direction for a new substituent on an atom that
    % already has EXISTINGDIRS ([Kx3] unit vectors toward its current
    % neighbours, K = 0, 1, 2, 3, or more). Not a hybridization solver (no
    % attempt to tell sp/sp2/sp3 apart) -- this always targets a
    % tetrahedral (109.47 degree) spread, which is "wrong but harmless" for
    % a genuinely sp/sp2 centre (refine the angle afterward, e.g. by
    % dragging in the 3D view, or via "Optimize structure") but, more
    % importantly, is never degenerate: an earlier, simpler version of
    % this ("negate the sum of existing directions", special-casing only
    % the very first addition to a 109.5 degree tilt) kept every new
    % direction inside the same plane spanned by whatever came before,
    % which is fine for 2 substituents but leaves no third dimension to
    % place a 3rd/4th one into -- confirmed the hard way: adding 4
    % hydrogens one at a time to the same carbon left 3 of them exactly
    % coplanar, and the geometric mean of two of those three coincided
    % with the third, so two hydrogens landed on the exact same point.
    %
    % This version instead finds the rigid rotation (Kabsch's algorithm,
    % via SVD -- the same method used to best-fit one point set onto
    % another) that maps a fixed, genuinely 3D reference tetrahedron's
    % first K vertex directions onto EXISTINGDIRS as closely as possible,
    % then returns that same rotation applied to the reference
    % tetrahedron's (K+1)-th vertex. For K=1 this still has one free
    % rotational degree about the existing bond's own axis (not pinned
    % down by a single point correspondence) -- any one valid completion
    % is fine, since what matters is that it is never coplanar with
    % whatever gets added next, which is guaranteed here because the
    % reference tetrahedron's own 4 vertices never are. Verified
    % directly: chaining K=0,1,2,3 calls on a fresh atom reproduces a
    % perfect tetrahedron (all 6 pairwise angles exactly 109.47 degrees).
    %
    % A 5th+ substituent (K>=4, chemically unusual already) has no
    % tetrahedral slot left to align to; it falls back to pointing away
    % from the mean of all existing directions instead.
        idealTetrahedron = [1 1 1; 1 -1 -1; -1 1 -1; -1 -1 1] / sqrt(3);
        K = size(existingDirs, 1);
        if K == 0
            dirVec = idealTetrahedron(1,:);
            return
        end
        if K >= 4
            s = sum(existingDirs, 1);
            if norm(s) < 1e-6
                dirVec = cross(existingDirs(1,:), existingDirs(2,:));
            else
                dirVec = -s;
            end
            dirVec = dirVec / norm(dirVec);
            return
        end
        H = idealTetrahedron(1:K,:)' * existingDirs;
        [U, ~, V] = svd(H);
        detSign = sign(det(V * U'));
        R = V * diag([1 1 detSign]) * U';
        dirVec = (R * idealTetrahedron(K+1,:)')';
        dirVec = dirVec / norm(dirVec);
    end

% -------------------------------------------------------------------------
    function R = rotationAligningVectors(u, w)
    % Returns the rotation matrix R (3x3, so that for a row vector p,
    % P_ROTATED = p*R') with R*u(:) == w(:), for unit row vectors U, W.
    % Closed form (cross-product/Rodrigues variant, avoids normalizing the
    % rotation axis explicitly): with cr = cross(u,w), c = dot(u,w), and
    % K = skew(cr), R = I + K + K^2*(1-c)/||cr||^2 -- equivalent to the
    % standard sin/cos Rodrigues formula since ||cr|| = sin(theta) and
    % c = cos(theta) for unit U, W.
        c = dot(u, w);
        cr = cross(u, w);
        s2 = dot(cr, cr);
        if s2 < 1e-12
            if c > 0
                R = eye(3);
            else
                % U and W point in exactly opposite directions: rotate 180
                % degrees about any axis perpendicular to U.
                perp = cross(u, [1 0 0]);
                if norm(perp) < 1e-6
                    perp = cross(u, [0 1 0]);
                end
                n = perp / norm(perp);
                R = 2 * (n' * n) - eye(3);
            end
            return
        end
        K = [0 -cr(3) cr(2); cr(3) 0 -cr(1); -cr(2) cr(1) 0];
        R = eye(3) + K + K*K*((1 - c) / s2);
    end

% -------------------------------------------------------------------------
    function onRemoveAtom()
        if selAtomIdx < 1
            return
        end
        idx = selAtomIdx;
        symbols(idx) = [];
        xyz(idx, :) = [];
        if ~isempty(bondList)
            mask = bondList(:,1) == idx | bondList(:,2) == idx;
            bondList(mask, :) = [];
            bondList(bondList(:,1) > idx, 1) = bondList(bondList(:,1) > idx, 1) - 1;
            bondList(bondList(:,2) > idx, 2) = bondList(bondList(:,2) > idx, 2) - 1;
        end
        selAtomIdx = 0;
        lastAtomsSelection = [];
        updateAtomsTable();
        updateBondsTable();
        redraw();
        updateButtonStates();
        snapshotPush();
        statusLabel.Text = sprintf('Removed atom %d.', idx);
    end

% -------------------------------------------------------------------------
    function onBondPair()
        if numel(lastAtomsSelection) ~= 2
            return
        end
        a = lastAtomsSelection(1);
        b = lastAtomsSelection(2);
        pairMask = (bondList(:,1)==a & bondList(:,2)==b) | (bondList(:,1)==b & bondList(:,2)==a);
        if any(pairMask)
            bondList(pairMask, :) = [];
            statusLabel.Text = sprintf('Removed bond between atoms %d and %d.', a, b);
        else
            bondList(end+1, :) = [a, b, 1];
            statusLabel.Text = sprintf('Added bond between atoms %d and %d.', a, b);
        end
        updateBondsTable();
        redraw();
        snapshotPush();
    end

% -------------------------------------------------------------------------
    function xyzOut = rigidAlignPairToVector(xyzIn, aIdx, bIdx, targetDir)
    % Rotates+translates XYZIN (every atom, rigidly -- one rotation
    % matrix applied to the whole set) so that atom BIDX - atom AIDX
    % points exactly along TARGETDIR (any nonzero 3-vector; not
    % necessarily unit length, and not necessarily a coordinate axis --
    % that is the one difference from an earlier, axis-only version of
    % this, which special-cased X/Y/Z; this is the same algorithm,
    % generalized). Atom AIDX itself ends up exactly ON the line through
    % the origin in direction TARGETDIR (its component perpendicular to
    % TARGETDIR is zeroed); BIDX lands on that same line too "for free",
    % since after the rotation AIDX->BIDX is an exact multiple of
    % TARGETDIR, so BIDX's perpendicular component already equals
    % AIDX's before this translation. Any other atom that was already
    % collinear with AIDX/BIDX beforehand lands on the line too, for the
    % same reason (collinearity with a fixed line is itself invariant
    % under a rigid motion) -- same "for free" logic as
    % ONALIGNTOPLANE's coplanar case. Never changes any bond length or
    % angle. Throws (caught by the caller, shown via UIALERT) if the two
    % atoms coincide or TARGETDIR is the zero vector.
        v = xyzIn(bIdx,:) - xyzIn(aIdx,:);
        vLen = norm(v);
        if vLen < 1e-9
            error('G_molecule_builder:coincidentAtoms', ...
                'The two selected atoms are at the same position; there is no direction to align.');
        end
        dLen = norm(targetDir);
        if dLen < 1e-9
            error('G_molecule_builder:zeroVector', ...
                'The target vector is zero; there is no direction to align to.');
        end
        dHat = targetDir / dLen;
        R = rotationAligningVectors(v / vLen, dHat);
        pivot = xyzIn(aIdx,:);
        xyzRot = (xyzIn - pivot) * R' + pivot;
        perp = xyzRot(aIdx,:) - dot(xyzRot(aIdx,:), dHat) * dHat;
        xyzOut = xyzRot - perp;
    end

% -------------------------------------------------------------------------
    function onAlignToAxis()
    % A1 = lastAtomsSelection(1) (lower row index -- UITABLE multi-selection
    % is always reported sorted, so "A1"/"A2" here is an index convention,
    % not click order; same convention ONBONDPAIR already uses). See
    % RIGIDALIGNPAIRTOVECTOR for the algorithm (a coordinate axis is just
    % one particular TARGETDIR).
        if numel(lastAtomsSelection) ~= 2
            return
        end
        a = lastAtomsSelection(1);
        b = lastAtomsSelection(2);
        switch alignAxisDD.Value
            case 'X', targetAxis = [1 0 0];
            case 'Y', targetAxis = [0 1 0];
            otherwise, targetAxis = [0 0 1];
        end
        try
            xyz = rigidAlignPairToVector(xyz, a, b, targetAxis);
        catch ME
            uialert(fig, ME.message, 'Cannot align');
            return
        end
        showAxesCheck.Value = true;
        updateAtomsTable();
        redraw();
        snapshotPush();
        statusLabel.Text = sprintf('Aligned atoms %d->%d onto the %s axis.', a, b, alignAxisDD.Value);
    end

% -------------------------------------------------------------------------
    function onAlignToVector()
    % Same as ONALIGNTOAXIS, but the target direction is an arbitrary
    % vector typed into the Vec x/y/z fields instead of a coordinate
    % axis -- e.g. the molecule's own dipole moment vector (as reported
    % by a Gaussian 'opt'/'freq' job on it), to view/orient the molecule
    % relative to that direction. See RIGIDALIGNPAIRTOVECTOR for the
    % algorithm.
        if numel(lastAtomsSelection) ~= 2
            return
        end
        a = lastAtomsSelection(1);
        b = lastAtomsSelection(2);
        targetDir = [vecXField.Value, vecYField.Value, vecZField.Value];
        try
            xyz = rigidAlignPairToVector(xyz, a, b, targetDir);
        catch ME
            uialert(fig, ME.message, 'Cannot align');
            return
        end
        showAxesCheck.Value = true;
        updateAtomsTable();
        redraw();
        snapshotPush();
        statusLabel.Text = sprintf('Aligned atoms %d->%d onto vector [%.4g %.4g %.4g].', ...
            a, b, targetDir(1), targetDir(2), targetDir(3));
    end

% -------------------------------------------------------------------------
    function onAlignToPlane()
    % A1/A2/A3 = lastAtomsSelection, sorted ascending (same row-index
    % convention as ONALIGNTOAXIS/ONBONDPAIR). The plane through the three
    % selected atoms has normal CROSS(P2-P1, P3-P1); rotating that normal
    % onto the chosen coordinate plane's own normal (Z for XY, Y for XZ,
    % X for YZ) makes the atoms' plane parallel to the target one, then a
    % single-coordinate translation (analogous to ONALIGNTOAXIS's two-
    % coordinate one) zeroes that one coordinate for A1 -- A2 and A3 land
    % exactly on it too "for free", since all three remain coplanar
    % throughout and that plane is now exactly the target coordinate
    % plane. Only ever a rigid rotation+translation of the WHOLE molecule
    % (no projection, no change to any bond length/angle): if the
    % molecule (or some larger subset of it) was already coplanar with
    % A1/A2/A3 before this -- e.g. a deliberately built planar fragment --
    % every atom in that coplanar subset also lands exactly on the target
    % plane, not just A1/A2/A3 themselves. An atom that was NOT already
    % coplanar with A1/A2/A3 is carried along rigidly and generally will
    % not end up on the plane -- picking A1/A2/A3 so that the atoms you
    % actually want flattened are coplanar with them is up to the user.
        if numel(lastAtomsSelection) ~= 3
            return
        end
        a = lastAtomsSelection(1);
        b = lastAtomsSelection(2);
        c = lastAtomsSelection(3);
        p1 = xyz(a,:); p2 = xyz(b,:); p3 = xyz(c,:);
        normalVec = cross(p2 - p1, p3 - p1);
        nLen = norm(normalVec);
        if nLen < 1e-9
            uialert(fig, ['The three selected atoms are collinear (or coincide); ', ...
                'they do not define a unique plane.'], 'Cannot align');
            return
        end
        switch alignPlaneDD.Value
            case 'XY', targetAxis = [0 0 1]; offComp = 3;
            case 'XZ', targetAxis = [0 1 0]; offComp = 2;
            otherwise, targetAxis = [1 0 0]; offComp = 1; % 'YZ'
        end
        R = rotationAligningVectors(normalVec / nLen, targetAxis);
        pivot = xyz(a,:);
        xyzRot = (xyz - pivot) * R' + pivot;
        shift = zeros(1, 3);
        shift(offComp) = -xyzRot(a, offComp);
        xyz = xyzRot + shift;
        showAxesCheck.Value = true;
        updateAtomsTable();
        redraw();
        snapshotPush();
        statusLabel.Text = sprintf('Aligned atoms %d,%d,%d onto the %s plane.', a, b, c, alignPlaneDD.Value);
    end

% -------------------------------------------------------------------------
    function onRemoveBond()
        if isempty(lastBondsSelection)
            return
        end
        bondList(lastBondsSelection, :) = [];
        lastBondsSelection = [];
        updateBondsTable();
        redraw();
        updateButtonStates();
        snapshotPush();
    end

% -------------------------------------------------------------------------
    function onAtomsTableEdit(evt)
    % The uitable widget itself already commits EVT.NEWDATA into
    % ATOMSTABLE.DATA before this callback runs (that's how a live,
    % interactive edit always behaves) -- but what it committed is the
    % raw typed value, which for the Symbol column can differ from the
    % normalized form actually stored in SYMBOLS (e.g. typed "cl" ->
    % normalized "Cl"). UPDATEATOMSTABLE() is still called below so the
    % displayed table always matches the real, normalized state exactly.
        row = evt.Indices(1,1);
        col = evt.Indices(1,2);
        switch col
            case 2
                newSym = normalizeSymbol(evt.NewData);
                if isempty(newSym)
                    atomsTable.Data{row,2} = evt.PreviousData;
                    return
                end
                symbols{row} = newSym;
            case 3
                xyz(row,1) = evt.NewData;
            case 4
                xyz(row,2) = evt.NewData;
            case 5
                xyz(row,3) = evt.NewData;
        end
        updateAtomsTable();
        redraw();
        snapshotPush();
    end

% -------------------------------------------------------------------------
    function onAtomsTableSelection(evt)
        if isempty(evt.Indices)
            lastAtomsSelection = [];
            selAtomIdx = 0;
        else
            rows = unique(evt.Indices(:,1));
            lastAtomsSelection = rows;
            if numel(rows) == 1
                selAtomIdx = rows(1);
            else
                selAtomIdx = 0;
            end
        end
        updateButtonStates();
        redraw();
    end

% -------------------------------------------------------------------------
    function onBondsTableEdit(evt)
    % See ONATOMSTABLEEDIT's comment: the widget already committed the
    % raw typed value before this callback runs, which can differ from
    % the rounded integer actually stored (e.g. typed "2.7" -> stored 3).
        row = evt.Indices(1,1);
        newOrder = round(evt.NewData);
        if ~ismember(newOrder, [1 2 3])
            bondsTable.Data{row,3} = evt.PreviousData;
            return
        end
        bondList(row,3) = newOrder;
        updateBondsTable();
        redraw();
        snapshotPush();
    end

% -------------------------------------------------------------------------
    function onBondsTableSelection(evt)
        if isempty(evt.Indices)
            lastBondsSelection = [];
        else
            lastBondsSelection = unique(evt.Indices(:,1));
        end
        updateButtonStates();
    end

% -------------------------------------------------------------------------
    function onLoadGjf()
        [f, p] = pickOpenFile({'*.gjf;*.com;*.in', 'Gaussian input (*.gjf, *.com, *.in)'}, ...
            'Load Gaussian input');
        if isequal(f, 0)
            return
        end
        fn = fullfile(p, f);
        try
            ginp = G_read_input(fn);
        catch ME
            uialert(fig, ME.message, 'G_read_input error');
            return
        end

        symbols = ginp.symbols(:);
        xyz = ginp.xyz;
        try
            bt = G16_get_bond_length(struct('symbols', {symbols}, 'xyz', xyz));
            bondList = [bt.Atom1, bt.Atom2, ones(height(bt), 1)];
        catch
            bondList = zeros(0, 3);
        end

        if ~isnan(ginp.charge), chargeField.Value = ginp.charge; end
        if ~isnan(ginp.mult),   multField.Value   = ginp.mult;   end
        if ~isempty(ginp.route), routeField.Value = ginp.route; end
        if ~isempty(ginp.title), titleField.Value = ginp.title; end
        if ~isempty(ginp.chk),   chkField.Value   = ginp.chk;   end
        if ~isempty(ginp.mem),   memField.Value   = ginp.mem;   end
        if ~isempty(ginp.nproc)
            v = str2double(ginp.nproc);
            if ~isnan(v), nprocField.Value = v; end
        end

        selAtomIdx = 0;
        lastAtomsSelection = [];
        lastBondsSelection = [];
        historyStack = {struct('symbols', {symbols}, 'xyz', xyz, 'bondList', bondList)};
        historyIdx = 1;

        updateAtomsTable();
        updateBondsTable();
        redraw();
        updateButtonStates();
        updateUndoRedoButtons();
        statusLabel.Text = sprintf('Loaded %s (%d atoms, %d bonds auto-detected).', fn, numel(symbols), size(bondList,1));
    end

% -------------------------------------------------------------------------
    function onOptimizeStructure()
    % Dispatches to one of two geometry clean-up methods, picked via
    % OPTMETHODDD: the built-in, dependency-free G_ROUGH_OPTIMIZE_GEOMETRY
    % (a simple geometric penalty -- not a real force field), or Open
    % Babel's own force-field minimizer (G_OPENBABEL_OPTIMIZE_GEOMETRY --
    % a real, validated force field, but requires Open Babel installed
    % separately; see that function's own error message if it is not).
    % Both take/return the same (symbols, xyz, bondList) shape, so adding
    % this dispatcher only needed one shared UI row instead of a growing
    % button per method.
        if numel(symbols) < 2
            uialert(fig, 'Add at least two bonded atoms first.', 'Nothing to optimize');
            return
        end
        method = optMethodDD.Value;
        statusLabel.Text = sprintf('Running %s...', method);
        drawnow;
        try
            switch method
                case 'Quick cleanup (built-in)'
                    xyz = G_rough_optimize_geometry(symbols, xyz, bondList);
                case 'Open Babel: MMFF94'
                    xyz = G_openbabel_optimize_geometry(symbols, xyz, bondList, 'ForceField', 'MMFF94');
                case 'Open Babel: MMFF94s'
                    xyz = G_openbabel_optimize_geometry(symbols, xyz, bondList, 'ForceField', 'MMFF94s');
                case 'Open Babel: UFF'
                    xyz = G_openbabel_optimize_geometry(symbols, xyz, bondList, 'ForceField', 'UFF');
            end
        catch ME
            uialert(fig, ME.message, sprintf('%s error', method));
            statusLabel.Text = sprintf('%s failed (see error).', method);
            return
        end
        updateAtomsTable();
        redraw();
        snapshotPush();
        statusLabel.Text = sprintf('%s applied. Refine further by hand, or run a real "opt" job in Gaussian.', method);
    end

% -------------------------------------------------------------------------
    function onGenerateInput()
        if isempty(symbols)
            uialert(fig, 'Add at least one atom first.', 'Nothing to export');
            return
        end
        [f, p] = pickSaveFile({'*.gjf', 'Gaussian input (*.gjf)'}, 'Save Gaussian input as', 'molecule.gjf');
        if isequal(f, 0)
            return
        end
        outFile = fullfile(p, f);

        nvArgs = {'Title', titleField.Value};
        if ~isempty(strtrim(chkField.Value)), nvArgs = [nvArgs, {'Chk', strtrim(chkField.Value)}]; end
        if ~isempty(strtrim(memField.Value)), nvArgs = [nvArgs, {'Mem', strtrim(memField.Value)}]; end
        if nprocField.Value > 0,              nvArgs = [nvArgs, {'Nproc', nprocField.Value}];      end

        try
            G_write_gjf(outFile, symbols, xyz, routeField.Value, chargeField.Value, multField.Value, nvArgs{:});
            statusLabel.Text = sprintf('Gaussian input saved: %s', outFile);
        catch ME
            uialert(fig, ME.message, 'G_write_gjf error');
        end
    end

end % function G_molecule_builder
