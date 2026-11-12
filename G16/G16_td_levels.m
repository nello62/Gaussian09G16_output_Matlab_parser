function L = G16_td_levels(logfile, varargin)
%G16_TD_LEVELS  Orbital-level diagram of a single TD-DFT calculation, with
%   arrows for the dominant excitations of selected excited states.
%
%   L = G16_td_levels('NMA_wB97XD_augTZ_TD_trans.log')
%   L = G16_td_levels('NMA_WB97XD_DEF2TZVP_TD_ALL.LOG', 'Job', 1, ...
%                     'States', [3 6 8], 'Labels', {19,'\pi_{nb}'; 20,'n_O'; 22,'\pi^*'})
%
%   The diagram shows the Kohn-Sham orbital energies (eV) around the
%   frontier: occupied levels dark grey, virtual levels light grey, thicker
%   for HOMO, LUMO and the orbitals involved in the drawn excitations.
%   Energies are drawn at their true values; near-degenerate levels are not
%   shifted, only their labels are spread apart.
%   For every selected state one arrow is drawn per dominant excitation
%   i -> a (weight 2c^2 >= 'MinWeight'): one colour per state, line width
%   growing with the weight, weight in % written along the arrow and the
%   state name (S1, S6, ...) at the arrow tip. The legend gives energy,
%   wavelength and oscillator strength of each state.
%
%   Name-value options
%     'Job'       : job of a --Link1-- chain (default: last job with excited states)
%     'States'    : 'auto' (default) = first singlet + the 'NBright' most
%                   intense singlets up to 'EMax'; or a vector of state numbers.
%                   Same rule and colours as G16_td_plot: the states and
%                   colours of the two figures match.
%     'NBright'   : number of intense states chosen by 'auto' (default 2)
%     'EMax'      : upper energy (eV) for 'auto' (default: first singlet + 2.5 eV)
%     'MinWeight' : smallest weight drawn as an arrow (default 0.15)
%     'NOcc'      : occupied orbitals shown, HOMO included (default: enough to
%                   include all drawn excitations, at least 3)
%     'NVirt'     : virtual orbitals shown, LUMO included (default: as NOcc)
%     'Labels'    : chemical labels {MO, 'label'; ...}, TeX allowed
%     'Title'     : figure title (default: file name and job)
%     'YLim'      : energy range [Emin Emax] in eV (default: automatic)
%     'Plot'      : true (default) | false (data only)
%
%   Output L (struct): file, job, HOMO, occ_eV, virt_eV (all alpha orbital
%   energies, eV), states (state numbers drawn), arrows (struct array with
%   state, mult, from, to, weight, E_eV, nm, f).
%
%   Requires G16_td_plot.m (reads the excited states).

p = inputParser;
p.addRequired('logfile', @(s) ischar(s) || isstring(s));
p.addParameter('Job', []);
p.addParameter('States', 'auto');
p.addParameter('NBright', 2);
p.addParameter('EMax', []);
p.addParameter('MinWeight', 0.15);
p.addParameter('NOcc', []);
p.addParameter('NVirt', []);
p.addParameter('Labels', {});
p.addParameter('Title', '');
p.addParameter('YLim', []);
p.addParameter('Plot', true);
p.parse(logfile, varargin{:});
o = p.Results;
logfile = char(logfile);

% ---------------- excited states and orbital energies ----------------
args = {'Plot', false, 'MinWeight', min(o.MinWeight, 0.05), 'States', o.States, ...
        'NBright', o.NBright, 'EMax', o.EMax};
if ~isempty(o.Job), args = [args, {'Job', o.Job}]; end
T = G16_td_plot(logfile, args{:});          % same state selection as the spectrum
[occ, virt] = read_orbital_energies(logfile, T.job);
H = T.HOMO;  if isnan(H), H = numel(occ); end
E = [occ(:); virt(:)];

% ---------------- states to draw ----------------
pick = find(ismember(T.state, T.selected));
if isnumeric(o.States)
    miss = setdiff(o.States, T.state);
    if ~isempty(miss), warning('G16_td_levels: states %s not found.', mat2str(miss)); end
end
ar = struct('state', {}, 'mult', {}, 'from', {}, 'to', {}, 'weight', {}, 'E_eV', {}, 'nm', {}, 'f', {});
for k = pick
    tr = T.transitions{k};
    for m = 1:numel(tr)
        if tr(m).deexc || tr(m).weight < o.MinWeight, continue, end
        ar(end+1) = struct('state', T.state(k), 'mult', T.mult{k}, 'from', tr(m).from, 'to', tr(m).to, ...
            'weight', tr(m).weight, 'E_eV', T.E_eV(k), 'nm', T.nm(k), 'f', T.f(k)); %#ok<AGROW>
    end
end
L = struct('file', logfile, 'job', T.job, 'HOMO', H, 'occ_eV', occ, 'virt_eV', virt, ...
           'states', T.state(pick), 'arrows', ar);

% printout
fprintf('\n-- G16_td_levels: HOMO = %d (%.2f eV), LUMO = %d (%.2f eV) --\n', H, E(H), H + 1, E(H + 1));
for k = 1:numel(ar)
    fprintf('  %s%-3d %6.3f eV  f = %.4f   %3d -> %-3d  (%s -> %s)  %3.0f%%\n', ...
        state_prefix(ar(k).mult), ar(k).state, ar(k).E_eV, ar(k).f, ar(k).from, ar(k).to, ...
        hl_label(ar(k).from, H), hl_label(ar(k).to, H), 100 * ar(k).weight);
end
if ~o.Plot, return; end

% ---------------- orbitals shown ----------------
involved = unique([[ar.from], [ar.to]]);
nOcc = o.NOcc;  nVirt = o.NVirt;
if isempty(nOcc),  nOcc  = max(3, H - min([involved(involved <= H), H]) + 1); end
if isempty(nVirt), nVirt = max(3, max([involved(involved > H), H + 1]) - H); end
mo  = max(1, H - nOcc + 1):min(numel(E), H + nVirt);
Emo = E(mo);

% ---------------- plot ----------------
pal = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; 0.600 0.318 0.784; ...
       0.839 0.208 0.420; 0.000 0.600 0.700; 0.700 0.550 0.100; 0.400 0.400 0.400];
[~, fname, fext] = fileparts(logfile);
ttl = char(o.Title);  if isempty(ttl), ttl = sprintf('%s%s, job %d', fname, fext, T.job); end
figure('Color', 'w', 'Name', ttl);
ax = axes;  hold(ax, 'on');
W = 6;                                      % width of the level lines (axis units)

% levels
for k = 1:numel(mo)
    isOcc = mo(k) <= H;
    isImp = ismember(mo(k), involved) || mo(k) == H || mo(k) == H + 1;
    col = [0.25 0.25 0.25];  if ~isOcc, col = [0.55 0.55 0.55]; end
    lw  = 0.8;  if isImp, lw = 2.4; end
    plot(ax, [0 W], [Emo(k) Emo(k)], '-', 'Color', col, 'LineWidth', lw);
end

% labels of the important levels, spread apart if too close
yr  = [min(Emo) max(Emo)];
sel = find(arrayfun(@(m) ismember(m, involved) || m == H || m == H + 1, mo));
[ys, ord] = sort(Emo(sel));  sel = sel(ord);
yl  = spread_labels(ys(:)', 0.035 * diff(yr));
for k = 1:numel(sel)
    m = mo(sel(k));
    s = sprintf('%d (%s)', m, hl_label(m, H));
    cl = chem_label(o.Labels, m);
    if ~isempty(cl), s = [s '  ' cl]; end %#ok<AGROW>
    text(ax, W + 0.25, yl(k), s, 'FontSize', 9, 'VerticalAlignment', 'middle', 'Color', [0.15 0.15 0.15]);
    if abs(yl(k) - ys(k)) > 1e-6
        plot(ax, [W, W + 0.22], [ys(k) yl(k)], ':', 'Color', [0.6 0.6 0.6]);
    end
end
text(ax, -0.2, mean(occ(max(1, H - nOcc + 1):H)), 'occupied', 'Rotation', 90, ...
     'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
text(ax, -0.2, mean(E(H + 1:mo(end))), 'virtual', 'Rotation', 90, ...
     'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);

% arrows, grouped by state
na = numel(ar);
hLeg = gobjects(0);  legTxt = {};
if na > 0
    stl = T.state(pick);                     % colour index = position in the selection
    xs  = W * (1:na) / (na + 1);
    for k = 1:na
        ci  = find(stl == ar(k).state, 1);
        col = pal(mod(ci - 1, size(pal, 1)) + 1, :);
        y1 = E(ar(k).from);  y2 = E(ar(k).to);
        quiver(ax, xs(k), y1, 0, y2 - y1, 0, 'Color', col, 'LineWidth', 0.8 + 2.4 * ar(k).weight, ...
               'MaxHeadSize', 0.5 / max(abs(y2 - y1), 1));
        text(ax, xs(k), y1 + 0.45 * (y2 - y1), sprintf('%.0f%%', 100 * ar(k).weight), 'Rotation', 90, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', col);
        text(ax, xs(k), y2 + 0.012 * diff(yr), sprintf('%s%d', state_prefix(ar(k).mult), ar(k).state), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, ...
             'FontWeight', 'bold', 'Color', col);
    end
    for s = stl
        if ~any([ar.state] == s), continue, end
        ci  = find(stl == s, 1);
        col = pal(mod(ci - 1, size(pal, 1)) + 1, :);
        a1  = ar(find([ar.state] == s, 1));
        hLeg(end+1) = plot(ax, nan, nan, '-', 'Color', col, 'LineWidth', 2.5); %#ok<AGROW>
        legTxt{end+1} = sprintf('%s%d   %.2f eV (%.0f nm)   f = %.4f', ...
            state_prefix(a1.mult), s, a1.E_eV, a1.nm, a1.f); %#ok<AGROW>
    end
end

% axes
if isempty(o.YLim)
    pad = 0.06 * diff(yr);  ylim(ax, [yr(1) - pad, yr(2) + 2 * pad]);
else
    ylim(ax, o.YLim);
end
xlim(ax, [-0.6, W + 3.2]);
set(ax, 'XTick', [], 'XColor', 'none', 'FontSize', 10, 'Box', 'off');
ylabel(ax, 'Orbital energy \epsilon (eV)');
title(ax, ttl, 'Interpreter', 'none', 'FontWeight', 'normal');
if ~isempty(hLeg)
    legend(ax, hLeg, legTxt, 'Location', 'southoutside', 'Box', 'off', 'FontSize', 9);
end
end

% ===========================================================================
function [occ, virt] = read_orbital_energies(file, job)
% Alpha orbital energies (eV) from the ground-state (SCF density) population
% analysis of the given job.
txt   = fileread(file);
lines = regexp(txt, '\r?\n', 'split');
starts = find(~cellfun(@isempty, strfind(lines, 'Entering Link 1')));
if isempty(starts), starts = 1; end
starts(1) = 1;
ends = [starts(2:end) - 1, numel(lines)];
blk = lines(starts(job):ends(job));
k0 = find(~cellfun(@isempty, strfind(blk, 'Population analysis using the SCF')), 1, 'last');
if isempty(k0), k0 = 1; end
occ = [];  virt = [];  started = false;
for k = k0:numel(blk)
    s = blk{k};
    if ~isempty(strfind(s, 'Alpha  occ. eigenvalues'))
        occ = [occ, parse_vals(s)];  started = true; %#ok<AGROW>
    elseif ~isempty(strfind(s, 'Alpha virt. eigenvalues'))
        virt = [virt, parse_vals(s)];  started = true; %#ok<AGROW>
    elseif started
        break
    end
end
if isempty(occ)
    error('G16_td_levels: orbital energies not found in job %d of %s.', job, file);
end
occ = occ * 27.211386;  virt = virt * 27.211386;
end

function v = parse_vals(s)
% values after '--' (Gaussian may glue negative numbers, e.g. -10.1234-10.0012)
s = s(strfind(s, '--') + 2:end);
v = str2double(regexp(s, '-?\d+\.\d+', 'match'));
end

function y = spread_labels(y, dmin)
% push label positions apart (minimum distance dmin), keeping their order
for it = 1:100
    moved = false;
    for k = 2:numel(y)
        if y(k) - y(k-1) < dmin
            mid = (y(k) + y(k-1)) / 2;
            y(k-1) = mid - dmin / 2;  y(k) = mid + dmin / 2;
            moved = true;
        end
    end
    if ~moved, break, end
end
end

function s = hl_label(m, H)
if m <= H
    d = H - m;  if d == 0, s = 'H'; else, s = sprintf('H-%d', d); end
else
    d = m - H - 1;  if d == 0, s = 'L'; else, s = sprintf('L+%d', d); end
end
end

function s = state_prefix(mult)
switch lower(mult)
    case 'singlet', s = 'S';
    case 'triplet', s = 'T';
    otherwise,      s = '';
end
end

function s = chem_label(lab, m)
s = '';
if isempty(lab), return, end
for k = 1:size(lab, 1)
    if isequal(lab{k, 1}, m), s = lab{k, 2}; return, end
end
end
