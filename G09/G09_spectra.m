function sp = G09_spectra(filename, varargin)
% G09_SPECTRA  Extracts IR and Raman spectra from a Gaussian 09 output file
%              and generates Lorentzian-broadened continuous spectra.
%
%   sp = G09_SPECTRA(filename)
%   sp = G09_SPECTRA(filename, Name, Value, ...)
%
%   Difference from G16: Gaussian 09 writes only ONE Harmonic frequencies
%   section (even for opt+freq jobs), so no 'section' parameter is needed.
%
%   Optional parameters (Name-Value):
%       'FWHM'      - Lorentzian full width at half maximum (cm-1, default: 10)
%       'xmin'      - lower wavenumber limit (cm-1, default: 0)
%       'xmax'      - upper wavenumber limit (cm-1, default: 4000)
%       'dx'        - grid step (cm-1, default: 1)
%       'normalize' - normalise continua to maximum = 1 (default: false)
%       'plot'      - generate figure (default: false)
%       'Lines'     - pre-read cell array of file lines (from
%                     G09_READ_LINES), to skip re-reading the file when it
%                     has already been read elsewhere (e.g. G09_READ_ALL).
%                     Default {} (read the file normally).
%
%   OUTPUT  struct sp with fields:
%       .freq        [Nmodes x 1]   frequencies (cm-1)
%       .IR          [Nmodes x 1]   IR intensities (KM/Mole)
%       .Raman       [Nmodes x 1]   Raman activities (A^4/AMU), [] if absent
%       .has_Raman   logical
%       .Nmodes      int
%       .x           [Ngrid x 1]    wavenumber grid (cm-1)
%       .IR_cont     [Ngrid x 1]    Lorentzian-broadened IR spectrum
%       .Raman_cont  [Ngrid x 1]    Lorentzian-broadened Raman spectrum
%       .FWHM        double
%       .filename    char

% -------------------------------------------------------------------------
% Parse arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename',  @ischar);
addParameter(p, 'FWHM',      10,    @isnumeric);
addParameter(p, 'xmin',      0,     @isnumeric);
addParameter(p, 'xmax',      4000,  @isnumeric);
addParameter(p, 'dx',        1,     @isnumeric);
addParameter(p, 'normalize', false, @islogical);
addParameter(p, 'plot',      false, @islogical);
addParameter(p, 'Lines',     {},    @iscell);
parse(p, filename, varargin{:});

FWHM    = p.Results.FWHM;
xmin    = p.Results.xmin;
xmax    = p.Results.xmax;
dx      = p.Results.dx;
do_norm = p.Results.normalize;
do_plot = p.Results.plot;

% -------------------------------------------------------------------------
% Read file
% -------------------------------------------------------------------------
lines = p.Results.Lines;
if isempty(lines)
    lines = G09_read_lines(filename);
end
N     = numel(lines);

% -------------------------------------------------------------------------
% Find the single "Harmonic frequencies" section
% -------------------------------------------------------------------------
sec_starts = find(~cellfun(@isempty, strfind(lines, 'Harmonic frequencies')));

if isempty(sec_starts)
    error('G09_spectra: no "Harmonic frequencies" section found in %s', filename);
end
sec_start = sec_starts(end);   % use last (only one in G09, but safe)

% -------------------------------------------------------------------------
% Parse Frequencies, IR Inten, Raman Activ
% -------------------------------------------------------------------------
freqs  = [];
IRs    = [];
Ramans = [];

end_pat = regexp_end_section();

for k = sec_start : N
    ln = lines{k};
    if k > sec_start && ~isempty(regexp(ln, end_pat, 'once')), break; end

    if ~isempty(regexp(ln, '^\s*Frequencies\s*--', 'once'))
        freqs = [freqs; sscanf(parse_rhs(ln), '%f')]; %#ok<AGROW>
        continue
    end
    if ~isempty(regexp(ln, '^\s*IR Inten\s*--', 'once'))
        IRs = [IRs; sscanf(parse_rhs(ln), '%f')]; %#ok<AGROW>
        continue
    end
    if ~isempty(regexp(ln, '^\s*Raman Activ\s*--', 'once'))
        Ramans = [Ramans; sscanf(parse_rhs(ln), '%f')]; %#ok<AGROW>
        continue
    end
end

if isempty(freqs)
    error('G09_spectra: no frequencies read from %s', filename);
end

Nmodes    = numel(freqs);
has_Raman = numel(Ramans) == Nmodes;

fix_vec = @(v) [v(1:min(end,Nmodes)); zeros(max(0,Nmodes-numel(v)),1)];
IRs = fix_vec(IRs);
if ~has_Raman, Ramans = []; end

% -------------------------------------------------------------------------
% Lorentzian-broadened continua
% -------------------------------------------------------------------------
x     = (xmin : dx : xmax)';
gamma = FWHM / 2;

IR_cont    = zeros(size(x));
Raman_cont = zeros(size(x));

for m = 1:Nmodes
    L = gamma^2 ./ ((x - freqs(m)).^2 + gamma^2);
    IR_cont = IR_cont + IRs(m) * L;
    if has_Raman
        Raman_cont = Raman_cont + Ramans(m) * L;
    end
end

if do_norm
    if max(IR_cont)    > 0, IR_cont    = IR_cont    / max(IR_cont);    end
    if max(Raman_cont) > 0, Raman_cont = Raman_cont / max(Raman_cont); end
end

if ~has_Raman, Raman_cont = []; end

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
sp.freq        = freqs;
sp.IR          = IRs;
sp.Raman       = Ramans;
sp.has_Raman   = has_Raman;
sp.Nmodes      = Nmodes;
sp.x           = x;
sp.IR_cont     = IR_cont;
sp.Raman_cont  = Raman_cont;
sp.FWHM        = FWHM;
sp.filename    = filename;

% -------------------------------------------------------------------------
% Optional plot
% -------------------------------------------------------------------------
if do_plot
    plot_spectra_internal(sp);
end

end  % G09_spectra


% =========================================================================
%  Local functions
% =========================================================================
function s = parse_rhs(ln)
idx = strfind(ln, '--');
s   = ln(idx(1)+2 : end);
end

function pat = regexp_end_section()
pat = '^\s*(-{20,}|Thermochemistry|Zero-point|Normal termination|Leave Link)';
end

function plot_spectra_internal(sp)
[~, fname] = fileparts(sp.filename);
nrows = 1 + sp.has_Raman;
fig   = figure('Color', 'white', 'Name', fname, 'NumberTitle', 'off');

if sp.has_Raman
    ax1 = subplot(2, 1, 1, 'Parent', fig);
    hold(ax1, 'on');
    for m = 1:sp.Nmodes
        if sp.Raman(m) > 0
            line(ax1, [sp.freq(m) sp.freq(m)], [0 sp.Raman(m)], ...
                 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
    end
    plot(ax1, sp.x, sp.Raman_cont, 'Color', [0.15 0.45 0.80], 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Raman  FWHM=%g cm^{-1}', sp.FWHM));
    set(ax1, 'XDir', 'reverse', 'Box', 'on');
    xlabel(ax1, 'Wavenumber (cm^{-1})', 'FontSize', 10);
    ylabel(ax1, 'Raman activity (A^4 AMU^{-1})', 'FontSize', 10);
    title(ax1, [strrep(fname,'_','\_'), '  —  Raman'], 'FontSize', 10);
    legend(ax1, 'show', 'Location', 'northeast', 'Box', 'off');
    xlim(ax1, [sp.x(1) sp.x(end)]);
    ax_ir = subplot(2, 1, 2, 'Parent', fig);
else
    ax_ir = axes('Parent', fig);
end

hold(ax_ir, 'on');
for m = 1:sp.Nmodes
    if sp.IR(m) > 0
        line(ax_ir, [sp.freq(m) sp.freq(m)], [0 sp.IR(m)], ...
             'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end
plot(ax_ir, sp.x, sp.IR_cont, 'Color', [0.85 0.20 0.15], 'LineWidth', 1.5, ...
     'DisplayName', sprintf('IR  FWHM=%g cm^{-1}', sp.FWHM));
set(ax_ir, 'XDir', 'reverse', 'Box', 'on');
xlabel(ax_ir, 'Wavenumber (cm^{-1})', 'FontSize', 10);
ylabel(ax_ir, 'IR intensity (KM mol^{-1})', 'FontSize', 10);
title(ax_ir, [strrep(fname,'_','\_'), '  —  IR'], 'FontSize', 10);
legend(ax_ir, 'show', 'Location', 'northeast', 'Box', 'off');
xlim(ax_ir, [sp.x(1) sp.x(end)]);
end
