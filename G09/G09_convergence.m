function cv = G09_convergence(filename, varargin)
% G09_CONVERGENCE  Extracts geometry optimisation convergence criteria
%                  from a Gaussian 09 output file.
%
%   cv = G09_CONVERGENCE(filename)
%   cv = G09_CONVERGENCE(filename, 'plot', true)
%
%   The convergence block format is identical to G16.
%
%   OUTPUT  struct cv with fields:
%       .MaxForce       [Nsteps x 1]   maximum force (a.u.)
%       .RMSForce       [Nsteps x 1]   RMS force (a.u.)
%       .MaxDisp        [Nsteps x 1]   maximum displacement (a.u.)
%       .RMSDisp        [Nsteps x 1]   RMS displacement (a.u.)
%       .thr_MaxForce   double
%       .thr_RMSForce   double
%       .thr_MaxDisp    double
%       .thr_RMSDisp    double
%       .converged      logical
%       .conv_step      int or NaN
%       .Nsteps         int
%       .filename       char

p = inputParser;
addRequired(p,  'filename', @ischar);
addParameter(p, 'plot',     false, @islogical);
parse(p, filename, varargin{:});
do_plot = p.Results.plot;

lines = G09_read_lines(filename);
N     = numel(lines);

pat = '([\d.]+(?:D[+-]?\d+)?)\s+([\d.]+(?:D[+-]?\d+)?)\s+(YES|NO)';

MaxForce = []; RMSForce = []; MaxDisp  = []; RMSDisp  = [];
thrs = [NaN NaN NaN NaN];
converged  = false; conv_step  = NaN; step_count = 0;

k = 1;
while k <= N
    ln = lines{k};

    if ~isempty(strfind(ln, 'Item               Value     Threshold  Converged?'))
        if k+4 > N, break; end
        t1 = parse_conv_line(lines{k+1}, pat);
        t2 = parse_conv_line(lines{k+2}, pat);
        t3 = parse_conv_line(lines{k+3}, pat);
        t4 = parse_conv_line(lines{k+4}, pat);

        if ~isempty(t1) && ~isempty(t2) && ~isempty(t3) && ~isempty(t4)
            step_count = step_count + 1;
            MaxForce(end+1) = t1(1); %#ok<AGROW>
            RMSForce(end+1) = t2(1); %#ok<AGROW>
            MaxDisp(end+1)  = t3(1); %#ok<AGROW>
            RMSDisp(end+1)  = t4(1); %#ok<AGROW>
            if isnan(thrs(1))
                thrs = [t1(2) t2(2) t3(2) t4(2)];
            end
            all_yes = t1(3) && t2(3) && t3(3) && t4(3);
            if all_yes && ~converged
                converged = true;
                conv_step = step_count;
            end
        end
        k = k + 5; continue
    end
    k = k + 1;
end

if step_count == 0
    error('G09_convergence: no convergence blocks found in %s', filename);
end

cv.MaxForce     = MaxForce(:);
cv.RMSForce     = RMSForce(:);
cv.MaxDisp      = MaxDisp(:);
cv.RMSDisp      = RMSDisp(:);
cv.thr_MaxForce = thrs(1); cv.thr_RMSForce = thrs(2);
cv.thr_MaxDisp  = thrs(3); cv.thr_RMSDisp  = thrs(4);
cv.converged    = converged;
cv.conv_step    = conv_step;
cv.Nsteps       = step_count;
cv.filename     = filename;

fprintf('\nG09_convergence: %s\n', filename);
fprintf('  Steps read : %d\n', step_count);
if converged
    fprintf('  Converged  : YES  (step %d)\n', conv_step);
else
    fprintf('  Converged  : NO\n');
end
fprintf('\n');

if do_plot && step_count > 0
    plot_convergence_internal(cv);
end

end  % G09_convergence


function out = parse_conv_line(ln, pat)
tok = regexp(ln, pat, 'tokens', 'once');
if isempty(tok)
    out = [];
    return
end
v1  = str2double(strrep(tok{1}, 'D', 'e'));
v2  = str2double(strrep(tok{2}, 'D', 'e'));
yes = strcmp(tok{3}, 'YES');
out = [v1, v2, double(yes)];
end


function plot_convergence_internal(cv)
[~, fname] = fileparts(cv.filename);
steps  = (1:cv.Nsteps)';
labels = {'Max Force','RMS Force','Max Displacement','RMS Displacement'};
data   = [cv.MaxForce, cv.RMSForce, cv.MaxDisp, cv.RMSDisp];
thrs   = [cv.thr_MaxForce, cv.thr_RMSForce, cv.thr_MaxDisp, cv.thr_RMSDisp];
colors = [0.15 0.45 0.80; 0.85 0.20 0.15; 0.10 0.65 0.30; 0.90 0.55 0.00];

fig = figure('Color','white','Name',fname,'NumberTitle','off');
for pi = 1:4
    ax = subplot(2,2,pi);
    hold(ax,'on');
    semilogy(ax, steps, data(:,pi), 'o-', 'Color', colors(pi,:), ...
             'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', colors(pi,:));
    yline(ax, thrs(pi), '--k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    if cv.converged && cv.conv_step <= cv.Nsteps
        semilogy(ax, cv.conv_step, data(cv.conv_step,pi), 'p', 'MarkerSize', 10, ...
                 'MarkerFaceColor', [0.2 0.8 0.2], 'MarkerEdgeColor', 'k', ...
                 'HandleVisibility', 'off');
    end
    set(ax,'Box','on','YGrid','on','XGrid','off');
    xlabel(ax,'Opt step','FontSize',9); ylabel(ax,labels{pi},'FontSize',9);
    title(ax,labels{pi},'FontSize',10); xlim(ax,[1 max(steps)]);
end
sgtitle(strrep(fname,'_','\_'),'FontSize',11,'Interpreter','tex');
end
