function outfile = G09_animate_mode(mol, nm, mode_idx, varargin)
% G09_ANIMATE_MODE  Exports an MP4 animation of a vibrational mode.
%
%   G09_ANIMATE_MODE(mol, nm, mode_idx) oscillates the molecule along the
%   mode's displacement vector (equilibrium +/- amplitude, like GaussView's
%   mode animations) and saves the result as an MP4 video (VideoWriter,
%   MPEG-4 profile).
%
%   OUTFILE = G09_ANIMATE_MODE(...) returns the path written to.
%
%   Name-Value parameters:
%       'Filename'        - output path (default: '<source>_mode<N>.mp4';
%                            '.mp4' is appended if missing)
%       'Scale'           - displacement amplitude scale, same meaning as
%                            G09_draw_mode's 'Scale' (default 1.5)
%       'FlipSign'        - invert the displacement direction (default false)
%       'AtomScale'       - see G09_draw_molecule (default 0.35)
%       'BondTol'         - see G09_draw_molecule (default 1.30)
%       'ShowLabels'      - see G09_draw_molecule (default false)
%       'FramesPerCycle'  - frames per oscillation period (default 30)
%       'NCycles'         - number of periods rendered (default 2)
%       'FPS'             - video frame rate (default 20)
%       'View'            - [azimuth elevation] in degrees, the starting
%                            camera orientation (default [] = MATLAB's
%                            standard 3D default view, azimuth=-37.5,
%                            elevation=30). Pass the output of MATLAB's
%                            own view(ax) to match a figure you have
%                            already rotated interactively.
%       'ShowWaitbar'     - print rendering progress to the command window
%                           (default true). A graphical waitbar is
%                           deliberately not used: MATLAB's waitbar()
%                           itself was found to error intermittently
%                           ("Not enough input arguments" inside its own
%                           title() call) when updated across many
%                           iterations on some systems/versions.
%
%   A figure window opens and animates live while the video is rendered;
%   avoid interacting with it until G09_ANIMATE_MODE returns.
%
%   Example:
%       mol = G09_structure('V_E00t.out');
%       nm  = G09_nmodes('V_E00t.out');
%       G09_animate_mode(mol, nm, 70, 'Scale', 2, 'FPS', 25);
%
%   See also G09_DRAW_MODE, G09_DRAW_MOLECULE, G09_MODEVIEWER.

% -------------------------------------------------------------------------
% Parse input arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'mol',       @isstruct);
addRequired(p,  'nm',        @isstruct);
addRequired(p,  'mode_idx',  @isnumeric);
addParameter(p, 'Filename',       '',    @ischar);
addParameter(p, 'Scale',          1.5,   @isnumeric);
addParameter(p, 'FlipSign',       false, @islogical);
addParameter(p, 'AtomScale',      0.35,  @isnumeric);
addParameter(p, 'BondTol',        1.30,  @isnumeric);
addParameter(p, 'ShowLabels',     false, @islogical);
addParameter(p, 'FramesPerCycle', 30,    @isnumeric);
addParameter(p, 'NCycles',        2,     @isnumeric);
addParameter(p, 'FPS',            20,    @isnumeric);
addParameter(p, 'View',           [],    @isnumeric);
addParameter(p, 'ShowWaitbar',    true,  @islogical);
parse(p, mol, nm, mode_idx, varargin{:});

scale            = p.Results.Scale;
flip_sign        = p.Results.FlipSign;
atom_scale       = p.Results.AtomScale;
bond_tol         = p.Results.BondTol;
show_labels      = p.Results.ShowLabels;
view_angle       = p.Results.View;
frames_per_cycle = round(p.Results.FramesPerCycle);
n_cycles         = round(p.Results.NCycles);
fps              = p.Results.FPS;
show_wb          = p.Results.ShowWaitbar;
outfile          = p.Results.Filename;

if mode_idx < 1 || mode_idx > nm.Nmodes
    error('G09_animate_mode:badIndex', ...
        'mode index %d is out of range [1, %d].', mode_idx, nm.Nmodes);
end
if mol.Natoms ~= nm.Natoms
    error('G09_animate_mode:sizeMismatch', ...
        'mol.Natoms (%d) does not match nm.Natoms (%d).', mol.Natoms, nm.Natoms);
end

% -------------------------------------------------------------------------
% Output filename
% -------------------------------------------------------------------------
if isempty(outfile)
    if isfield(mol, 'filename') && ~isempty(mol.filename)
        [~, fn] = fileparts(mol.filename);
    else
        fn = 'molecule';
    end
    outfile = sprintf('%s_mode%d.mp4', fn, mode_idx);
end
[~, ~, ext] = fileparts(outfile);
if ~strcmpi(ext, '.mp4')
    outfile = [outfile '.mp4'];
end

% -------------------------------------------------------------------------
% Displacement amplitude (mirrors G09_draw_mode's normalisation)
% -------------------------------------------------------------------------
U = nm.disp(:, :, mode_idx);
if flip_sign
    U = -U;
end
norms_i  = sqrt(sum(U.^2, 2));
max_norm = max(norms_i);
if max_norm == 0
    error('G09_animate_mode:zeroDisplacement', ...
        'Zero displacement vectors for mode %d.', mode_idx);
end
U_scaled = U / max_norm * scale;

% -------------------------------------------------------------------------
% Fixed axis limits across the whole oscillation, so the camera/box does
% not jitter frame to frame (G09_draw_molecule's own 'axis tight' would
% otherwise re-fit to each frame's instantaneous, slightly different extent)
% -------------------------------------------------------------------------
pad = 1.0;   % Angstrom, roughly one atom radius, so spheres are not clipped
extreme_pts = [mol.xyz - abs(U_scaled); mol.xyz + abs(U_scaled)];
xlim_fixed = [min(extreme_pts(:,1))-pad, max(extreme_pts(:,1))+pad];
ylim_fixed = [min(extreme_pts(:,2))-pad, max(extreme_pts(:,2))+pad];
zlim_fixed = [min(extreme_pts(:,3))-pad, max(extreme_pts(:,3))+pad];

% -------------------------------------------------------------------------
% Fixed bond list from the equilibrium geometry, so bonds do not
% appear/disappear frame to frame as instantaneous distances oscillate
% across the BondTol threshold (G09_draw_molecule's default distance-based
% detection would otherwise re-evaluate connectivity on every frame).
% -------------------------------------------------------------------------
bondTable = G09_get_bond_length(mol, 'Tolerance', bond_tol, 'IncludeH', true);
bond_order = zeros(height(bondTable), 1);
for bb = 1:height(bondTable)
    bond_order(bb) = local_classify_bond_order(char(bondTable.Sym1(bb)), ...
        char(bondTable.Sym2(bb)), bondTable.Distance_Ang(bb));
end
bond_list = [bondTable.Atom1, bondTable.Atom2, bond_order];

freq_str = sprintf('Mode %d - %.1f cm^{-1}', mode_idx, nm.freq(mode_idx));
if isfield(mol, 'filename') && ~isempty(mol.filename)
    [~, fn] = fileparts(mol.filename);
    fig_title = sprintf('%s\n%s', strrep(fn, '_', '\_'), freq_str);
else
    fig_title = freq_str;
end

% -------------------------------------------------------------------------
% Render frames and write the video
% -------------------------------------------------------------------------
fig = figure('Color', 'white', 'Name', sprintf('Animating Mode %d...', mode_idx), ...
             'NumberTitle', 'off');
ax = axes('Parent', fig);

vw = VideoWriter(outfile, 'MPEG-4');
vw.FrameRate = fps;
vw.Quality   = 90;
open(vw);

total_frames = frames_per_cycle * n_cycles;

if show_wb
    fprintf('G09_animate_mode: rendering mode %d animation (%d frames)...\n', ...
        mode_idx, total_frames);
end
progress_step = max(1, round(total_frames / 10));   % ~10 console updates total

mol_frame = mol;
cleaner = onCleanup(@() close(fig));   % ensures the figure closes even on error

for k = 0:total_frames-1
    phase = sin(2*pi*k/frames_per_cycle);
    mol_frame.xyz = mol.xyz + phase * U_scaled;

    % cla(ax) was found to leave stale Surface/Line/Text objects behind
    % in this 3D+lighting configuration (confirmed via findall: counts
    % kept growing frame to frame even immediately after cla), causing
    % the video to show all frames superimposed instead of one at a
    % time. Deleting and recreating the axes is the reliable fix.
    delete(ax);
    ax = axes('Parent', fig);
    G09_draw_molecule(mol_frame, 'Ax', ax, 'AtomScale', atom_scale, ...
        'BondTol', bond_tol, 'ShowLabels', show_labels, ...
        'ShowLegend', false, 'Title', fig_title, 'BondList', bond_list);
    if ~isempty(view_angle)
        view(ax, view_angle(1), view_angle(2));   % overrides G09_draw_molecule's default view(3)
    end
    xlim(ax, xlim_fixed);
    ylim(ax, ylim_fixed);
    zlim(ax, zlim_fixed);

    drawnow;
    frame = getframe(fig);
    writeVideo(vw, frame);

    if show_wb && (mod(k+1, progress_step) == 0 || k+1 == total_frames)
        fprintf('  %3.0f%%  (%d/%d frames)\n', 100*(k+1)/total_frames, k+1, total_frames);
    end
end
close(vw);

fprintf('G09_animate_mode: animation saved to %s (%d frames, %.1f fps)\n', ...
    outfile, total_frames, fps);

end % G09_animate_mode


% =========================================================================
function order = local_classify_bond_order(sym_i, sym_j, d)
%LOCAL_CLASSIFY_BOND_ORDER  Same C-C/C-N/C-O bond-length classification as
%   G09_draw_molecule's local classify_bond_order, duplicated here so the
%   bond order can be fixed once from the equilibrium geometry (see above)
%   instead of being re-evaluated from the oscillating instantaneous
%   distance on every animation frame.
    p = sort({upper(sym_i), upper(sym_j)});
    pair = [p{1}, p{2}];
    switch pair
        case 'CC'
            thresh = [1.27, 1.36];    % [triple/double, double/single]; the
                                      % double/single boundary is set below
                                      % the ~1.39-1.40 A aromatic C-C range,
                                      % so symmetric aromatic rings render
                                      % as all-single, not all-double.
        case 'CN'
            thresh = [1.22, 1.375];
        case 'CO'
            thresh = [1.165, 1.315];
        otherwise
            order = 1;
            return
    end
    if d < thresh(1)
        order = 3;
    elseif d < thresh(2)
        order = 2;
    else
        order = 1;
    end
end
