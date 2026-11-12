function T = G16_td_plot(logfile, varargin)
%G16_TD_PLOT  Reads the TD-DFT (or CIS) excited states of a Gaussian 16 log and
%   plots them in order of energy with their oscillator strengths.
%
%   T = G16_td_plot('NMA_wB97XD_def2TZVP_TD_trans.log')
%   T = G16_td_plot('NMA_WB97XD_DEF2TZVP_TD_ALL.LOG', 'Job', 2)
%   T = G16_td_plot(file, 'XUnit', 'nm', 'FWHM', 0.4, 'Label', 5)
%
%   Figure
%     top    : stick spectrum, oscillator strength f versus excitation
%              energy; singlets (or all states of an unrestricted run) as
%              stems, triplets (f = 0 by symmetry) as triangles on the axis.
%              The selected states ('States') are drawn in colour and
%              labelled with state number, energy and dominant orbital
%              excitation (H = HOMO, L = LUMO, e.g. "S6  H-1->L+1 77%").
%     bottom : (only if FWHM > 0) Gaussian-broadened absorption spectrum
%              epsilon (L mol^-1 cm^-1), with the same x axis.
%
%   Name-value options
%     'Job'       : job of a --Link1-- chain to read (1-based); default:
%                   the last job that contains excited states
%     'XUnit'     : 'eV' (default) | 'nm' | 'cm-1'
%     'FWHM'      : full width at half maximum of the Gaussian broadening,
%                   in eV (default 0.3; 0 = sticks only)
%     'States'    : states highlighted and labelled: 'auto' (default) = first
%                   singlet + the 'NBright' most intense singlets up to 'EMax';
%                   or a vector of state numbers. Same rule and same colours
%                   as G16_td_levels, so the two figures match.
%     'NBright'   : number of intense states chosen by 'auto' (default 2)
%     'EMax'      : upper energy (eV) for 'auto' (default: first singlet + 2.5 eV)
%     'NLabel'    : (old option) label the NLabel most intense states instead
%     'Triplets'  : true (default) | false, show the triplet states
%     'MinWeight' : smallest weight 2*c^2 of an orbital excitation kept in
%                   T.transitions (default 0.05)
%     'Overlay'   : false (default) | true: one single plot, with the
%                   broadened spectrum drawn as a shaded envelope behind the
%                   sticks, scaled so that its maximum equals the largest f
%                   (shape only; use Overlay = false for epsilon in units)
%     'Plot'      : true (default) | false
%     'Title'     : figure title (default: file name, job)
%
%   Output T (struct)
%     state, mult ('Singlet','Triplet','Doublet',...), sym, E_eV, nm,
%     cm1, f, S2       per-state arrays (in Gaussian's order = energy order)
%     transitions      cell array, per state: struct array with fields
%                      from, to, coef, weight (= 2 c^2 for closed shell),
%                      deexc (true for "<-" de-excitation lines), label
%     main             per-state string with the dominant excitation
%     HOMO             index of the HOMO (number of alpha electrons)
%     selected         state numbers highlighted (in order of energy)
%     job, filename
%
%   The excited-state block read is the LAST one in the job (for TD
%   optimisations or state-specific solvation this is the final one).
%
%   Example
%     T = G16_td_plot('NMA_wB97XD_def2TZVP_TD_trans.log','FWHM',0.3);
%     tab = table(T.state', T.mult', T.E_eV', T.nm', T.f', T.main', ...
%                 'VariableNames', {'State','Mult','E_eV','nm','f','Main'})

p = inputParser;
p.addRequired('logfile', @(s) ischar(s) || isstring(s));
p.addParameter('Job', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('XUnit', 'eV', @(s) any(strcmpi(s, {'eV','nm','cm-1'})));
p.addParameter('FWHM', 0.3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('States', 'auto');
p.addParameter('NBright', 2, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('EMax', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('NLabel', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
p.addParameter('Triplets', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('MinWeight', 0.05, @(x) isnumeric(x) && isscalar(x));
p.addParameter('Plot', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Overlay', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Title', '', @(s) ischar(s) || isstring(s));
p.parse(logfile, varargin{:});
o = p.Results;
logfile = char(logfile);

txt   = fileread(logfile);
lines = regexp(txt, '\r?\n', 'split');

% ---------------- jobs of a --Link1-- chain ----------------
starts = find(~cellfun(@isempty, strfind(lines, 'Entering Link 1')));
if isempty(starts), starts = 1; end
starts(1) = 1;
ends = [starts(2:end) - 1, numel(lines)];
hasTD = false(1, numel(starts));
for j = 1:numel(starts)
    hasTD(j) = any(~cellfun(@isempty, regexp(lines(starts(j):ends(j)), '^\s*Excited State\s+\d+:', 'once')));
end
if ~any(hasTD)
    error('G16_td_plot: no "Excited State" lines found in %s.', logfile);
end
if isempty(o.Job)
    job = find(hasTD, 1, 'last');
else
    job = o.Job;
    if job < 1 || job > numel(starts), error('G16_td_plot: Job %d out of range [1, %d].', job, numel(starts)); end
    if ~hasTD(job), error('G16_td_plot: job %d contains no excited states (jobs with TD: %s).', job, mat2str(find(hasTD))); end
end
blk = lines(starts(job):ends(job));

% ---------------- HOMO = number of alpha electrons ----------------
HOMO = NaN;
tok = regexp(blk, '^\s*(\d+)\s+alpha electrons', 'tokens', 'once');
k = find(~cellfun(@isempty, tok), 1, 'last');
if ~isempty(k), HOMO = str2double(tok{k}{1}); end

% ---------------- last excited-state block ----------------
isES = find(~cellfun(@isempty, regexp(blk, '^\s*Excited State\s+1:', 'once')));
blk = blk(isES(end):end);
esRe = ['^\s*Excited State\s+(\d+):\s+(\S+)\s+(-?[\d.]+)\s+eV\s+(-?[\d.]+)\s+nm\s+' ...
        'f=\s*(-?[\d.]+)(?:\s+<S\*\*2>=\s*(-?[\d.]+))?'];
trRe = '^\s*(\d+)([AB]?)\s*(->|<-)\s*(\d+)([AB]?)\s+(-?[\d.]+)\s*$';

T = struct('state', [], 'mult', {{}}, 'sym', {{}}, 'E_eV', [], 'nm', [], 'cm1', [], ...
           'f', [], 'S2', [], 'transitions', {{}}, 'main', {{}}, 'HOMO', HOMO, ...
           'job', job, 'filename', logfile);
n = 0;  i = 1;
while i <= numel(blk)
    t = regexp(blk{i}, esRe, 'tokens', 'once');
    if isempty(t), i = i + 1; continue, end
    st = str2double(t{1});
    if st == 1 && n > 0, break; end      % safety: a further block would start here
    n = n + 1;
    ms = regexp(t{2}, '^([A-Za-z]+)-?(.*)$', 'tokens', 'once');
    T.state(n) = st;
    T.mult{n}  = ms{1};
    T.sym{n}   = ms{2};
    T.E_eV(n)  = str2double(t{3});
    T.nm(n)    = str2double(t{4});
    T.cm1(n)   = T.E_eV(n) * 8065.54429;
    T.f(n)     = str2double(t{5});
    if numel(t) >= 6 && ~isempty(t{6}), T.S2(n) = str2double(t{6}); else, T.S2(n) = NaN; end

    tr = struct('from', {}, 'to', {}, 'coef', {}, 'weight', {}, 'deexc', {}, 'label', {});
    i = i + 1;
    while i <= numel(blk)
        u = regexp(blk{i}, trRe, 'tokens', 'once');
        if isempty(u), break; end
        c  = str2double(u{6});
        w  = 2 * c^2;                     % closed shell: weight of the singlet-adapted pair
        if ~isempty(u{2}), w = c^2; end   % unrestricted (A/B spin labels)
        if w >= o.MinWeight
            m = numel(tr) + 1;
            tr(m).from   = str2double(u{1});
            tr(m).to     = str2double(u{4});
            tr(m).coef   = c;
            tr(m).weight = w;
            tr(m).deexc  = strcmp(u{3}, '<-');
            tr(m).label  = [orb_label(tr(m).from, HOMO, u{2}) '->' orb_label(tr(m).to, HOMO, u{5})];
        end
        i = i + 1;
    end
    T.transitions{n} = tr;
    if isempty(tr)
        T.main{n} = '';
    else
        [~, kk] = max([tr.weight]);
        T.main{n} = sprintf('%s %.0f%%', tr(kk).label, 100 * tr(kk).weight);
    end
end
if n == 0, error('G16_td_plot: excited states could not be parsed in %s.', logfile); end

% ---------------- states to highlight ----------------
isSg = ~strcmpi(T.mult, 'Triplet');
if ~isempty(o.NLabel)
    [fs, ord] = sort(T.f, 'descend');  ord = ord(fs > 0);
    pick = sort(ord(1:min(o.NLabel, numel(ord))));
elseif ischar(o.States) && strcmpi(o.States, 'auto')
    iS = find(isSg);
    if isempty(iS), iS = 1:n; end
    emax = o.EMax;  if isempty(emax), emax = T.E_eV(iS(1)) + 2.5; end
    iB = iS(T.E_eV(iS) <= emax);
    [~, ord] = sort(T.f(iB), 'descend');
    pick = sort(unique([iS(1), iB(ord(1:min(o.NBright, numel(ord))))]));
else
    pick = find(ismember(T.state, o.States));
end
T.selected = T.state(pick);

% ---------------- printout ----------------
[~, fname, fext] = fileparts(logfile);
fprintf('\n-- G16_td_plot: %s%s (job %d), %d states, HOMO = %d --\n', fname, fext, job, n, HOMO);
fprintf('%5s %-9s %8s %8s %8s  %s\n', 'State', 'Mult', 'E (eV)', 'nm', 'f', 'Main excitation');
for k = 1:n
    fprintf('%5d %-9s %8.4f %8.2f %8.4f  %s\n', T.state(k), T.mult{k}, T.E_eV(k), T.nm(k), T.f(k), T.main{k});
end

if ~o.Plot, return; end

% ---------------- plot ----------------
switch lower(o.XUnit)
    case 'ev',   x = T.E_eV; xl = 'Excitation energy (eV)';
    case 'nm',   x = T.nm;   xl = 'Wavelength (nm)';
    case 'cm-1', x = T.cm1;  xl = 'Wavenumber (cm^{-1})';
end
isT = strcmpi(T.mult, 'Triplet');
isS = ~isT;
cS = [0.45 0.45 0.45];       % singlets (grey; selected states in colour)
cT = [0.922 0.408 0.204];    % triplets
ttl = char(o.Title);
if isempty(ttl), ttl = sprintf('%s%s, job %d', fname, fext, job); end

figure('Color', 'w', 'Name', ttl);
overlay = logical(o.Overlay) && o.FWHM > 0;
if o.FWHM > 0 && ~overlay
    ax1 = subplot(2, 1, 1);
else
    ax1 = axes;
end
hold(ax1, 'on'); box(ax1, 'on');
h = gobjects(0); lg = {};
if overlay
    % broadened envelope in the background, scaled to the largest f
    [xg, eps_] = broadened(T, isS, o.FWHM, o.XUnit);
    fmx = max([T.f, 1e-3]);
    if max(eps_) > 0
        he = fill(ax1, [xg, fliplr(xg)], [eps_ / max(eps_) * fmx, zeros(size(xg))], ...
                  [0.81 0.88 0.96], 'EdgeColor', [0.55 0.70 0.90], 'LineWidth', 0.8);
        h(end+1) = he;  lg{end+1} = sprintf('Envelope, FWHM %.2f eV (scaled)', o.FWHM);
    end
end
if any(isS)
    hs = stem(ax1, x(isS), T.f(isS), 'filled', 'Color', cS, 'LineWidth', 1.6, 'MarkerSize', 4);
    h(end+1) = hs;  lg{end+1} = 'Singlets';
    if all(cellfun(@isempty, regexpi(T.mult(isS), 'singlet'))), lg{end} = 'Excited states'; end
end
if o.Triplets && any(isT)
    ht = plot(ax1, x(isT), zeros(1, nnz(isT)), 'v', 'Color', cT, 'MarkerFaceColor', cT, 'MarkerSize', 6);
    h(end+1) = ht;  lg{end+1} = 'Triplets (f = 0)';
end
fmax = max([T.f, 1e-3]);
ylim(ax1, [-0.04 * fmax, 1.25 * fmax]);
switch lower(o.XUnit)                      % margin: broadening + room for labels
    case 'ev',   mE = max(o.FWHM, 0.25);  xr = [min(T.E_eV) - mE, max(T.E_eV) + mE];
    case 'nm',   mE = max(o.FWHM, 0.25);  xr = sort(1239.84198 ./ [max(T.E_eV) + mE, max(min(T.E_eV) - mE, 0.5)]);
    case 'cm-1', mE = max(o.FWHM, 0.25);  xr = [min(T.E_eV) - mE, max(T.E_eV) + mE] * 8065.54429;
end
xlim(ax1, xr);
if strcmpi(o.XUnit, 'nm'), set(ax1, 'XDir', 'reverse'); end
ylabel(ax1, 'Oscillator strength f');
if o.FWHM <= 0 || overlay, xlabel(ax1, xl); end
title(ax1, ttl, 'Interpreter', 'none');
grid(ax1, 'on'); set(ax1, 'GridAlpha', 0.15, 'FontSize', 10);
if numel(h) > 1, legend(ax1, h, lg, 'Location', 'northwest', 'Box', 'off', 'AutoUpdate', 'off'); end

% selected states: coloured stems and labels (same colours as G16_td_levels)
pal = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; 0.600 0.318 0.784; ...
       0.839 0.208 0.420; 0.000 0.600 0.700; 0.700 0.550 0.100; 0.400 0.400 0.400];
plain = all(cellfun(@isempty, regexpi(T.mult, 'singlet|triplet')));
for j = 1:numel(pick)
    k = pick(j);
    col = pal(mod(j - 1, size(pal, 1)) + 1, :);
    pre = 'S'; if isT(k), pre = 'T'; end
    if plain, pre = ''; end
    plot(ax1, [x(k) x(k)], [0 T.f(k)], '-', 'Color', col, 'LineWidth', 2.6);
    plot(ax1, x(k), T.f(k), 'o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 5);
    s = sprintf('%s%d  %.2f eV\n%s', pre, T.state(k), T.E_eV(k), T.main{k});
    text(ax1, x(k), T.f(k) + 0.03 * fmax, s, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom', 'FontSize', 8, 'Interpreter', 'none', 'Color', col);
end

% broadened spectrum in a second panel
if o.FWHM > 0 && ~overlay
    ax2 = subplot(2, 1, 2);
    [xg, eps_] = broadened(T, isS, o.FWHM, o.XUnit);
    plot(ax2, xg, eps_, '-', 'Color', [0.165 0.471 0.839], 'LineWidth', 2);
    hold(ax2, 'on'); box(ax2, 'on');
    xlim(ax2, xlim(ax1));
    if strcmpi(o.XUnit, 'nm'), set(ax2, 'XDir', 'reverse'); end
    xlabel(ax2, xl);
    ylabel(ax2, '\epsilon (L mol^{-1} cm^{-1})');
    title(ax2, sprintf('Gaussian broadening, FWHM = %.2f eV', o.FWHM), 'FontWeight', 'normal');
    grid(ax2, 'on'); set(ax2, 'GridAlpha', 0.15, 'FontSize', 10);
    linkaxes([ax1 ax2], 'x');
end
end

% ---------------------------------------------------------------------------
function [xg, eps_] = broadened(T, isS, fwhm, xunit)
% Gaussian-broadened absorption spectrum (epsilon, L mol^-1 cm^-1) of the
% singlet states; each band normalised so that integral(eps d(nu~)) = f / 4.319e-9
sig  = fwhm / (2 * sqrt(2 * log(2)));                       % eV
Eg   = linspace(min(T.E_eV) - 4 * fwhm, max(T.E_eV) + 4 * fwhm, 2000);
sigc = sig * 8065.54429;                                     % cm-1
eps_ = zeros(size(Eg));
for k = find(isS & T.f > 0)
    eps_ = eps_ + f_to_eps(T.f(k), sigc) * exp(-(Eg - T.E_eV(k)).^2 / (2 * sig^2));
end
switch lower(xunit)
    case 'ev',   xg = Eg;
    case 'nm',   ok = Eg > 0.5;  Eg = Eg(ok);  eps_ = eps_(ok);  xg = 1239.84198 ./ Eg;
    case 'cm-1', xg = Eg * 8065.54429;
end
end

% ---------------------------------------------------------------------------
function s = orb_label(k, HOMO, spin)
% 'H-1', 'H', 'L', 'L+2' ... (plus spin label for unrestricted runs)
if isnan(HOMO)
    s = sprintf('%d', k);
elseif k <= HOMO
    d = HOMO - k;
    if d == 0, s = 'H'; else, s = sprintf('H-%d', d); end
else
    d = k - HOMO - 1;
    if d == 0, s = 'L'; else, s = sprintf('L+%d', d); end
end
if ~isempty(spin), s = [s '(' spin ')']; end
end

function a = f_to_eps(f, sigma_cm)
% Peak molar absorptivity of a Gaussian band with oscillator strength f and
% standard deviation sigma (cm-1):  integral(eps dnu) = f / 4.319e-9
a = f / 4.319e-9 / (sigma_cm * sqrt(2 * pi));
end
