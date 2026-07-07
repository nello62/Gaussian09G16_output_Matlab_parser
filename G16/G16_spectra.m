function sp = G16_spectra(filename, varargin)
% G16_SPECTRA  Extracts IR and Raman spectra from a Gaussian 16 .out/.log file
%              and generates Lorentzian-broadened continuous spectra.
%
%   sp = G16_SPECTRA(filename)
%   sp = G16_SPECTRA(filename, Name, Value, ...)
%
%   Optional parameters (Name-Value):
%       'FWHM'      - Lorentzian full width at half maximum (cm⁻¹)  (default: 10)
%       'xmin'      - lower x-axis limit in cm⁻¹               (default: 0)
%       'xmax'      - upper x-axis limit in cm⁻¹               (default: 4000)
%       'dx'        - grid step in cm⁻¹                          (default: 1)
%       'normalize' - normalise continua to maximum = 1            (default: false)
%       'plot'      - generate figure after extraction                 (default: false)
%       'section'   - 'last' use the last Harmonic freq section      (default: 'last')
%                     'first' use the first section
%
%   OUTPUT  struct sp with fields:
%       .freq        [Nmodes×1]   frequenze in cm⁻¹
%       .IR          [Nmodes×1]   IR intensities (KM/Mole)
%       .Raman       [Nmodes×1]   Raman activities (Å⁴/AMU), [] if absent
%       .Nmodes      int          number of normal modes
%       .has_Raman   logical
%       .x           [Ngrid×1]    griglia wavenumber (cm⁻¹)
%       .IR_cont     [Ngrid×1]    spettro IR continuum (Lorentziano)
%       .Raman_cont  [Ngrid×1]    spettro Raman continuum ([], se assente)
%       .FWHM        double       FWHM usata
%       .filename    char
%
%   Uso tipico:
%       sp = G16_spectra('zeatin.out');
%       plot(sp.x, sp.Raman_cont)
%
%       sp = G16_spectra('calc.out', 'FWHM', 15, 'xmin', 400, 'plot', true)
%
%       % used together with G16_structure in the same pipeline:
%       mol = G16_structure('calc.out');
%       sp  = G16_spectra('calc.out');

% -------------------------------------------------------------------------
% Parse input arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename',  @ischar);
addParameter(p, 'FWHM',      10,     @isnumeric);
addParameter(p, 'xmin',      0,      @isnumeric);
addParameter(p, 'xmax',      4000,   @isnumeric);
addParameter(p, 'dx',        1,      @isnumeric);
addParameter(p, 'normalize', false,  @islogical);
addParameter(p, 'plot',      false,  @islogical);
addParameter(p, 'section',   'last', @ischar);
addParameter(p, 'Lines',     {},     @iscell);
parse(p, filename, varargin{:});

FWHM      = p.Results.FWHM;
xmin      = p.Results.xmin;
xmax      = p.Results.xmax;
dx        = p.Results.dx;
do_norm   = p.Results.normalize;
do_plot   = p.Results.plot;
sec_req   = lower(p.Results.section);

% -------------------------------------------------------------------------
% Read file
% -------------------------------------------------------------------------
lines = p.Results.Lines;
if isempty(lines)
    if ~isfile(filename)
        error('G16_spectra: file not found: %s', filename);
    end
    fid  = fopen(filename, 'r');
    raw  = fread(fid, '*char')';
    fclose(fid);
    lines = strsplit(raw, newline);
end

% -------------------------------------------------------------------------
% Find sections "Harmonic frequencies ... Raman scattering" (with Raman)
% and  "Harmonic frequencies ... IR intensities" (IR only, no Raman)
% -------------------------------------------------------------------------
idx_raman_sec = find(~cellfun(@isempty, ...
    regexp(lines, 'Harmonic frequencies.*Raman scattering')));
idx_IR_sec    = find(~cellfun(@isempty, ...
    regexp(lines, 'Harmonic frequencies.*IR intensities.*KM')));
% merge and remove duplicates (wording varies across G16 revisions)
all_sec = sort(unique([idx_raman_sec, idx_IR_sec]));

if isempty(all_sec)
    error('G16_spectra: no "Harmonic frequencies" section found in %s', filename);
end

% Section selection
switch sec_req
    case 'last',  sec_start = all_sec(end);
    case 'first', sec_start = all_sec(1);
    otherwise,    error('G16_spectra: ''section'' deve essere ''last'' o ''first''.');
end

% -------------------------------------------------------------------------
% Parsing: raccoglie Frequencies, IR Inten, Raman Activ
% -------------------------------------------------------------------------
freqs  = [];
IRs    = [];
Ramans = [];

for k = sec_start : numel(lines)
    ln = lines{k};

    % End of section: long "---" line or end-of-block keyword
    if k > sec_start && ~isempty(regexp(ln, ...
            '^\s*(-{10,}|Thermochemistry|Zero-point|Normal termination)', 'once'))
        break
    end

    % Frequencies --   val1   val2   val3
    if ~isempty(regexp(ln, '^\s*Frequencies\s*--', 'once'))
        vals = sscanf(parse_after_dashdash(ln), '%f');
        freqs = [freqs; vals];  %#ok<AGROW>
        continue
    end

    % IR Inten    --   val1   val2   val3
    if ~isempty(regexp(ln, '^\s*IR Inten\s*--', 'once'))
        vals = sscanf(parse_after_dashdash(ln), '%f');
        IRs = [IRs; vals];  %#ok<AGROW>
        continue
    end

    % Raman Activ --   val1   val2   val3
    if ~isempty(regexp(ln, '^\s*Raman Activ\s*--', 'once'))
        vals = sscanf(parse_after_dashdash(ln), '%f');
        Ramans = [Ramans; vals];  %#ok<AGROW>
        continue
    end
end

% Cleanup: remove "--" read as NaN or spurious values
% sscanf on "-- val val val" skips "--" automatically (not numeric)
% but check that vectors are aligned
Nmodes = numel(freqs);
if Nmodes == 0
    error('G16_spectra: nessuna frequenza letta dalla sezione selezionata.');
end

% Allinea lunghezze (paranoia)
if numel(IRs) > Nmodes,    IRs    = IRs(1:Nmodes);    end
if numel(Ramans) > Nmodes, Ramans = Ramans(1:Nmodes); end
if numel(IRs) < Nmodes,    IRs    = [IRs; zeros(Nmodes-numel(IRs),1)]; end

has_Raman = numel(Ramans) == Nmodes && Nmodes > 0;

% -------------------------------------------------------------------------
% Grid and Lorentzian convolution
% -------------------------------------------------------------------------
x    = (xmin : dx : xmax)';
Ngrid = numel(x);
gamma = FWHM / 2;   % half-width at half-maximum

% Area-normalised Lorentzian: L(x,x0) = (gamma/pi) / ((x-x0)^2 + gamma^2)
% Here we use the peak-normalised version (not area-normalised):
%   L(x,x0) = gamma^2 / ((x-x0)^2 + gamma^2)
% so that each stick contributes its intensity at the peak

IR_cont    = zeros(Ngrid, 1);
Raman_cont = zeros(Ngrid, 1);

for m = 1 : Nmodes
    x0 = freqs(m);
    L  = (gamma^2) ./ ((x - x0).^2 + gamma^2);   % Lorentziana picco=1
    IR_cont = IR_cont + IRs(m) * L;
    if has_Raman
        Raman_cont = Raman_cont + Ramans(m) * L;
    end
end

% Optional normalisation
if do_norm
    if max(IR_cont) > 0,    IR_cont    = IR_cont    / max(IR_cont);    end
    if max(Raman_cont) > 0, Raman_cont = Raman_cont / max(Raman_cont); end
end

% -------------------------------------------------------------------------
% Build output struct
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

if ~has_Raman
    sp.Raman      = [];
    sp.Raman_cont = [];
end

% -------------------------------------------------------------------------
% Plot optional
% -------------------------------------------------------------------------
if do_plot
    G16_plot_spectra(sp);
end

end  % G16_spectra


% =========================================================================
%  Local function: extract the numeric part after '--'
% =========================================================================
function s = parse_after_dashdash(ln)
% Splits the line at the first '--' and returns everything after it.
% E.g.: ' Frequencies --     21.53   28.54' -> '     21.53   28.54'
idx = strfind(ln, '--');
if isempty(idx)
    s = ln;
else
    s = ln(idx(1)+2 : end);
end
end

% =========================================================================
%  Local function for plotting
% =========================================================================
function G16_plot_spectra(sp)

[~, fname] = fileparts(sp.filename);
fname_tex  = strrep(fname, '_', '\_');

if sp.has_Raman
    nrows = 2;
else
    nrows = 1;
end

fig = figure('Color', 'white', 'Name', fname, 'NumberTitle', 'off');

% ── Raman ─────────────────────────────────────────────────────────────────
if sp.has_Raman
    ax1 = subplot(nrows, 1, 1, 'Parent', fig);
    hold(ax1, 'on');

    % stick spectrum
    for m = 1:sp.Nmodes
        if sp.Raman(m) > 0
            line(ax1, [sp.freq(m) sp.freq(m)], [0 sp.Raman(m)], ...
                 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, ...
                 'HandleVisibility', 'off');
        end
    end

    % continuum
    plot(ax1, sp.x, sp.Raman_cont, 'Color', [0.15 0.45 0.80], 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Raman (FWHM = %g cm^{-1})', sp.FWHM));

    set(ax1, 'XDir', 'reverse', 'Box', 'on', 'XGrid', 'off', 'YGrid', 'off');
    xlabel(ax1, 'Wavenumber (cm^{-1})', 'FontSize', 10);
    ylabel(ax1, 'Raman activity (Å^4 AMU^{-1})', 'FontSize', 10);
    title(ax1, ['Raman — ' fname_tex], 'FontSize', 11, 'Interpreter', 'tex');
    legend(ax1, 'show', 'Location', 'northeast', 'Box', 'off');
    xlim(ax1, [sp.x(1) sp.x(end)]);
end

% ── IR ────────────────────────────────────────────────────────────────────
ax2 = subplot(nrows, 1, nrows, 'Parent', fig);
hold(ax2, 'on');

for m = 1:sp.Nmodes
    if sp.IR(m) > 0
        line(ax2, [sp.freq(m) sp.freq(m)], [0 sp.IR(m)], ...
             'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, ...
             'HandleVisibility', 'off');
    end
end

plot(ax2, sp.x, sp.IR_cont, 'Color', [0.85 0.20 0.15], 'LineWidth', 1.5, ...
     'DisplayName', sprintf('IR (FWHM = %g cm^{-1})', sp.FWHM));

set(ax2, 'XDir', 'reverse', 'Box', 'on');
xlabel(ax2, 'Wavenumber (cm^{-1})', 'FontSize', 10);
ylabel(ax2, 'IR intensity (KM mol^{-1})', 'FontSize', 10);
title(ax2, ['IR — ' fname_tex], 'FontSize', 11, 'Interpreter', 'tex');
legend(ax2, 'show', 'Location', 'northeast', 'Box', 'off');
xlim(ax2, [sp.x(1) sp.x(end)]);

end
