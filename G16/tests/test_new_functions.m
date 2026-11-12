function test_new_functions()
%TEST_NEW_FUNCTIONS  Smoke tests for G16_td_plot, G16_td_levels, split_link1_log
%   and G_relight (toolbox update of November 2026).
%
%   Run from the folder that contains this file, with the toolbox on the path:
%     addpath(genpath('<repo>/G16')); addpath(genpath('<repo>/G_Utility'));
%     test_new_functions
%
%   Test data (tests/data):
%     NMA_wB97XD_augTZ_TD_trans.log  single TD job, wB97X-D/aug-cc-pVTZ, 30 singlets
%     NMA_WB97XD_DEF2TZVP_TD_ALL.LOG two TD jobs chained with --Link1--
%                                    (20 singlets + 20 triplets each; %chk on
%                                    two lines because Gaussian wraps Link0 at 80 cols)

here = fileparts(mfilename('fullpath'));
d    = fullfile(here, 'data');
f1   = fullfile(d, 'NMA_wB97XD_augTZ_TD_trans.log');
f2   = fullfile(d, 'NMA_WB97XD_DEF2TZVP_TD_ALL.LOG');
nfail = 0;

% ---------------- G16_td_plot ----------------
T = G16_td_plot(f1, 'Plot', false);
nfail = nfail + check(numel(T.state) == 30,              'td_plot: 30 states');
nfail = nfail + check(abs(T.E_eV(1) - 5.9030) < 1e-4,    'td_plot: S1 = 5.9030 eV');
nfail = nfail + check(T.HOMO == 20,                      'td_plot: HOMO = 20');
nfail = nfail + check(isequal(T.selected(:)', [1 5 10]), 'td_plot: auto selection = [1 5 10]');
nfail = nfail + check(contains(T.main{1}, 'H->L+8'),     'td_plot: S1 main = H->L+8');

T2 = G16_td_plot(f2, 'Job', 2, 'Plot', false);
nfail = nfail + check(T2.job == 2,                       'td_plot: Job option on --Link1-- chain');
nfail = nfail + check(numel(T2.state) == 40,             'td_plot: 20 singlets + 20 triplets');
nfail = nfail + check(any(strcmp(T2.mult, 'Triplet')),   'td_plot: triplets parsed');

% plotting path (figure created, then closed)
T3 = G16_td_plot(f1, 'FWHM', 0.3, 'XUnit', 'nm'); close(gcf);
T4 = G16_td_plot(f1, 'Overlay', true);             close(gcf);
nfail = nfail + check(isequal(T3.selected, T4.selected), 'td_plot: plotting runs (nm, Overlay)');

% ---------------- G16_td_levels ----------------
L = G16_td_levels(f1, 'Plot', false);
nfail = nfail + check(L.HOMO == 20,                      'td_levels: HOMO = 20');
nfail = nfail + check(isequal(L.states(:)', [1 5 10]),   'td_levels: same states as td_plot');
nfail = nfail + check(any([L.arrows.to] == 29),          'td_levels: arrow to LUMO+8 (MO 29)');
L2 = G16_td_levels(f1, 'States', [1 3], 'Labels', {20, 'n_O'; 29, '\pi^*'}); close(gcf);
nfail = nfail + check(isequal(L2.states(:)', [1 3]),     'td_levels: explicit States + Labels');

% ---------------- split_link1_log ----------------
tmp = fullfile(tempdir, ['split_test_' char(java.util.UUID.randomUUID)]);
mkdir(tmp);
copyfile(f2, tmp);
files = split_link1_log(fullfile(tmp, 'NMA_WB97XD_DEF2TZVP_TD_ALL.LOG'), 'naming', 'chk');
nfail = nfail + check(numel(files) == 2,                 'split: 2 job logs');
[~, n1] = fileparts(files{1});
nfail = nfail + check(strcmp(n1, 'NMA_wB97XD_def2TZVP_TD_trans'), 'split: name from wrapped %chk');
Ts = G16_td_plot(files{2}, 'Plot', false);
nfail = nfail + check(numel(Ts.state) == 40,             'split: split log readable by G16_td_plot');
rmdir(tmp, 's');

% ---------------- G_relight ----------------
fig = figure('Visible', 'off'); ax = axes(fig);
surf(ax, peaks(20)); G_relight(ax, [-60 40]);
nfail = nfail + check(numel(findobj(ax, 'Type', 'light')) == 2, 'relight: two lights created');
G_relight(ax, [0 90]);
nfail = nfail + check(numel(findobj(ax, 'Type', 'light')) == 2, 'relight: lights re-used, not duplicated');
close(fig);

fprintf('\n%d check(s) failed.\n', nfail);
if nfail > 0, error('test_new_functions: %d failure(s).', nfail); end
end

function n = check(ok, msg)
if ok, fprintf('  PASS  %s\n', msg); n = 0;
else,  fprintf('  FAIL  %s\n', msg); n = 1; end
end
