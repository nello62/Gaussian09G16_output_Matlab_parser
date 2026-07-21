function cv = G16_convergence(filename, varargin)
% G16_CONVERGENCE  Extracts geometry optimisation convergence criteria
%                  from a Gaussian 16 .out/.log file.
%
%   cv = G16_CONVERGENCE(filename)
%   cv = G16_CONVERGENCE(filename, 'plot', true)
%
%   OUTPUT  struct cv with fields:
%       .MaxForce       [Nsteps x 1]   maximum force at each step (a.u.)
%       .RMSForce       [Nsteps x 1]   RMS force at each step (a.u.)
%       .MaxDisp        [Nsteps x 1]   maximum displacement at each step (a.u.)
%       .RMSDisp        [Nsteps x 1]   RMS displacement at each step (a.u.)
%       .thr_MaxForce   double         convergence threshold for maximum force
%       .thr_RMSForce   double         convergence threshold for RMS force
%       .thr_MaxDisp    double         convergence threshold for maximum displacement
%       .thr_RMSDisp    double         convergence threshold for RMS displacement
%       .converged      logical        true if all four criteria were satisfied
%       .conv_step      int            step at which convergence was reached (NaN if not)
%       .Nsteps         int            total number of convergence blocks read
%       .filename       char           source file path
%
%   Example:
%       cv = G16_convergence('V_E00t.out', 'plot', true)
%       cv.converged
%       cv.conv_step

% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename', @ischar);
addParameter(p, 'plot',     false, @islogical);
parse(p, filename, varargin{:});
do_plot = p.Results.plot;

% -------------------------------------------------------------------------
if ~isfile(filename)
    error('G16_convergence: file not found: %s', filename);
end
fid  = fopen(filename, 'r');
raw  = fread(fid, '*char')';
fclose(fid);
lines = strsplit(raw, newline);
G16_check_gaussian_match(lines, filename);
N = numel(lines);

% -------------------------------------------------------------------------
% Parsing
% -------------------------------------------------------------------------
% Block schema:
%         Item               Value     Threshold  Converged?
%  Maximum Force            0.001498     0.000450     NO
%  RMS     Force            0.000248     0.000300     YES
%  Maximum Displacement     0.094400     0.001800     NO
%  RMS     Displacement     0.020965     0.001200     NO

MaxForce = [];
RMSForce = [];
MaxDisp  = [];
RMSDisp  = [];

thr_MxF = NaN; thr_RMSF = NaN;
thr_MxD = NaN; thr_RMSD = NaN;

converged  = false;
conv_step  = NaN;
step_count = 0;

k = 1;
while k <= N
    ln = lines{k};

    % Detect convergence block header
    if ~isempty(strfind(ln, 'Item               Value     Threshold  Converged?'))
        % Read the 4 data lines that follow
        if k+4 > N, break; end

        ln1 = lines{k+1};  % Maximum Force
        ln2 = lines{k+2};  % RMS Force
        ln3 = lines{k+3};  % Maximum Displacement
        ln4 = lines{k+4};  % RMS Displacement

        t1 = parse_conv_line(ln1);
        t2 = parse_conv_line(ln2);
        t3 = parse_conv_line(ln3);
        t4 = parse_conv_line(ln4);

        if ~isempty(t1) && ~isempty(t2) && ~isempty(t3) && ~isempty(t4)
            step_count = step_count + 1;
            MaxForce(end+1) = t1(1); %#ok<AGROW>
            RMSForce(end+1) = t2(1); %#ok<AGROW>
            MaxDisp(end+1)  = t3(1); %#ok<AGROW>
            RMSDisp(end+1)  = t4(1); %#ok<AGROW>

            % Thresholds are constant; read them only on the first block
            if isnan(thr_MxF)
                thr_MxF  = t1(2);
                thr_RMSF = t2(2);
                thr_MxD  = t3(2);
                thr_RMSD = t4(2);
            end

            % Check if all four criteria are satisfied
            all_yes = t1(3) && t2(3) && t3(3) && t4(3);
            if all_yes && ~converged
                converged = true;
                conv_step = step_count;
            end
        end
        k = k + 5;
        continue
    end

    k = k + 1;
end

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
cv.MaxForce     = MaxForce(:);
cv.RMSForce     = RMSForce(:);
cv.MaxDisp      = MaxDisp(:);
cv.RMSDisp      = RMSDisp(:);
cv.thr_MaxForce = thr_MxF;
cv.thr_RMSForce = thr_RMSF;
cv.thr_MaxDisp  = thr_MxD;
cv.thr_RMSDisp  = thr_RMSD;
cv.converged    = converged;
cv.conv_step    = conv_step;
cv.Nsteps       = step_count;
cv.filename     = filename;

% -------------------------------------------------------------------------
% Print summary
% -------------------------------------------------------------------------
fprintf('\nG16_convergence: %s\n', filename);
fprintf('  Steps read  : %d\n', step_count);
if converged
    fprintf('  Converged   : YES  (step %d)\n', conv_step);
else
    fprintf('  Converged   : NO\n');
end
if step_count > 0
    fprintf('  Last step   : MaxF=%.2e (thr %.2e)  RMSF=%.2e (thr %.2e)\n', ...
        MaxForce(end), thr_MxF, RMSForce(end), thr_RMSF);
    fprintf('                MaxD=%.2e (thr %.2e)  RMSD=%.2e (thr %.2e)\n', ...
        MaxDisp(end),  thr_MxD, RMSDisp(end),  thr_RMSD);
end
fprintf('\n');

% -------------------------------------------------------------------------
% Plot optional
% -------------------------------------------------------------------------
if do_plot && step_count > 0
    plot_convergence(cv);
end

end  % G16_convergence


% =========================================================================
%  Local function: parse convergence line
%  Returns [value, threshold, converged_flag], or [] if line does not match
% =========================================================================
function out = parse_conv_line(ln)
% Format: " Maximum Force            0.001498     0.000450     NO "
%      o:  " RMS     Force            0.000248     0.000300     YES"
% Use regexp to capture both numeric values and the YES/NO flag

tok = regexp(ln, ...
    '([\d.]+(?:D[+-]?\d+)?)\s+([\d.]+(?:D[+-]?\d+)?)\s+(YES|NO)', ...
    'tokens', 'once');

if isempty(tok)
    out = [];
    return
end

% Convert Fortran D exponent notation (e.g. 1.234D-05) to MATLAB double
v1  = str2double(strrep(tok{1}, 'D', 'e'));
v2  = str2double(strrep(tok{2}, 'D', 'e'));
yes = strcmp(tok{3}, 'YES');

out = [v1, v2, double(yes)];
end


% =========================================================================
%  Convergence plot
% =========================================================================
function plot_convergence(cv)
[~, fname] = fileparts(cv.filename);
steps = (1:cv.Nsteps)';

colors = [0.15 0.45 0.80;   % blu   - MaxForce
          0.85 0.20 0.15;   % rosso - RMSForce
          0.10 0.65 0.30;   % verde - MaxDisp
          0.90 0.55 0.00];  % arancio - RMSDisp

labels    = {'Max Force', 'RMS Force', 'Max Displacement', 'RMS Displacement'};
data_all  = [cv.MaxForce, cv.RMSForce, cv.MaxDisp, cv.RMSDisp];
thrs      = [cv.thr_MaxForce, cv.thr_RMSForce, cv.thr_MaxDisp, cv.thr_RMSDisp];

fig = figure('Color', 'white', 'Name', fname, 'NumberTitle', 'off');

for pi = 1:4
    ax = subplot(2, 2, pi);
    hold(ax, 'on');

    semilogy(ax, steps, data_all(:, pi), ...
             'Color', colors(pi,:), 'LineWidth', 1.5, 'Marker', 'o', ...
             'MarkerSize', 4, 'MarkerFaceColor', colors(pi,:));

    % Soglia
    yline(ax, thrs(pi), '--k', 'LineWidth', 1.0, ...
          'HandleVisibility', 'off');

    % Convergence marker
    if cv.converged && cv.conv_step <= cv.Nsteps
        semilogy(ax, cv.conv_step, data_all(cv.conv_step, pi), ...
                 'p', 'MarkerSize', 10, 'MarkerFaceColor', [0.2 0.8 0.2], ...
                 'MarkerEdgeColor', [0 0 0], 'HandleVisibility', 'off');
    end

    set(ax, 'Box', 'on', 'YGrid', 'on', 'XGrid', 'off');
    xlabel(ax, 'Opt step', 'FontSize', 9);
    ylabel(ax, labels{pi}, 'FontSize', 9);
    title(ax,  labels{pi}, 'FontSize', 10);
    xlim(ax, [1, max(steps)]);
end

sgtitle(strrep(fname, '_', '\_'), 'FontSize', 11, 'Interpreter', 'tex');
end
