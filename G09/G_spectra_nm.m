function sp = G_spectra_nm(nm, varargin)
% G_SPECTRA_NM  Generates Lorentzian-broadened IR/Raman spectra from an
%               already-parsed nm (normal modes) struct, without reading
%               or re-parsing a Gaussian output file.
%
%   sp = G_SPECTRA_NM(nm)
%   sp = G_SPECTRA_NM(nm, Name, Value, ...)
%
%   Accepts the nm struct returned by G09_NMODES / G16_NMODES, or the
%   .nm sub-struct returned by G09_FCHK_READ / G16_FCHK_READ -- all four
%   share the same field layout (.freq, .IR, .Raman, .has_Raman,
%   .Nmodes, .filename), so this single function works with any of them
%   regardless of Gaussian version or data source (.out/.log vs .fchk).
%   In particular, this is what fills the gap left by G09_FCHK_READ /
%   G16_FCHK_READ: their output struct has a .nm field with frequencies
%   and IR/Raman intensities, but (unlike G09_SPECTRA/G16_SPECTRA) does
%   not itself build the broadened continuum spectrum -- this function
%   does exactly that, directly from data.nm, with no .out/.log file
%   needed at all.
%
%   Same broadening algorithm and Name-Value parameters as G09_SPECTRA /
%   G16_SPECTRA (peak-normalised Lorentzian convolution), just applied to
%   an already-parsed nm struct instead of reading a file:
%       'FWHM'      - Lorentzian full width at half maximum (cm^-1) (default: 10)
%       'xmin'      - lower x-axis limit in cm^-1               (default: 0)
%       'xmax'      - upper x-axis limit in cm^-1               (default: 4000)
%       'dx'        - grid step in cm^-1                        (default: 1)
%       'normalize' - normalise continua to maximum = 1         (default: false)
%       'plot'      - generate figure after extraction          (default: false)
%
%   OUTPUT  struct sp with fields:
%       .freq        [Nmodes x 1]  frequencies (cm^-1), copied from nm.freq
%       .IR          [Nmodes x 1]  IR intensities (KM/Mole)
%       .Raman       [Nmodes x 1]  Raman activities (Å^4/AMU), [] if absent
%       .Nmodes      int
%       .has_Raman   logical
%       .x           [Ngrid x 1]   wavenumber grid (cm^-1)
%       .IR_cont     [Ngrid x 1]   IR continuum spectrum (Lorentzian)
%       .Raman_cont  [Ngrid x 1]   Raman continuum spectrum ([] if absent)
%       .FWHM        double
%       .filename    char          copied from nm.filename if present, '' otherwise
%
%   Example:
%       data = G16_fchk_read('molecule.fchk');
%       sp   = G_spectra_nm(data.nm, 'FWHM', 15, 'plot', true);
%
%       nm = G16_nmodes('molecule.out');
%       sp = G_spectra_nm(nm, 'xmin', 400, 'normalize', true);
%
%   See also G09_SPECTRA, G16_SPECTRA, G09_FCHK_READ, G16_FCHK_READ,
%   G09_NMODES, G16_NMODES.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

% -------------------------------------------------------------------------
% Parse input arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'nm',        @isstruct);
addParameter(p, 'FWHM',      10,     @isnumeric);
addParameter(p, 'xmin',      0,      @isnumeric);
addParameter(p, 'xmax',      4000,   @isnumeric);
addParameter(p, 'dx',        1,      @isnumeric);
addParameter(p, 'normalize', false,  @islogical);
addParameter(p, 'plot',      false,  @islogical);
parse(p, nm, varargin{:});

FWHM    = p.Results.FWHM;
xmin    = p.Results.xmin;
xmax    = p.Results.xmax;
dx      = p.Results.dx;
do_norm = p.Results.normalize;
do_plot = p.Results.plot;

% -------------------------------------------------------------------------
% Validate input struct
% -------------------------------------------------------------------------
required = {'freq', 'IR', 'Nmodes'};
for k = 1:numel(required)
    if ~isfield(nm, required{k})
        error('G_spectra_nm: nm is missing field "%s". Use the nm struct returned by G09_nmodes/G16_nmodes, or the .nm field from G09_fchk_read/G16_fchk_read.', ...
            required{k});
    end
end

freqs  = nm.freq(:);
IRs    = nm.IR(:);
Nmodes = nm.Nmodes;

if Nmodes == 0 || isempty(freqs)
    error('G_spectra_nm: nm contains no vibrational modes (Nmodes = 0).');
end

has_Raman = isfield(nm, 'Raman') && ~isempty(nm.Raman) && numel(nm.Raman) == Nmodes;
if has_Raman
    Ramans = nm.Raman(:);
else
    Ramans = [];
end

if isfield(nm, 'filename')
    filename = nm.filename;
else
    filename = '';
end

% -------------------------------------------------------------------------
% Grid and Lorentzian convolution (identical algorithm to
% G09_SPECTRA/G16_SPECTRA: peak-normalised Lorentzian, not area-normalised)
% -------------------------------------------------------------------------
x     = (xmin : dx : xmax)';
Ngrid = numel(x);
gamma = FWHM / 2;   % half-width at half-maximum

IR_cont    = zeros(Ngrid, 1);
Raman_cont = zeros(Ngrid, 1);

for m = 1 : Nmodes
    x0 = freqs(m);
    L  = (gamma^2) ./ ((x - x0).^2 + gamma^2);   % Lorentzian, peak = 1
    IR_cont = IR_cont + IRs(m) * L;
    if has_Raman
        Raman_cont = Raman_cont + Ramans(m) * L;
    end
end

if do_norm
    if max(IR_cont) > 0, IR_cont = IR_cont / max(IR_cont); end
    if has_Raman && max(Raman_cont) > 0
        Raman_cont = Raman_cont / max(Raman_cont);
    end
end

% -------------------------------------------------------------------------
% Build output struct
% -------------------------------------------------------------------------
sp.freq       = freqs;
sp.IR         = IRs;
sp.Raman      = Ramans;
sp.has_Raman  = has_Raman;
sp.Nmodes     = Nmodes;
sp.x          = x;
sp.IR_cont    = IR_cont;
sp.Raman_cont = Raman_cont;
sp.FWHM       = FWHM;
sp.filename   = filename;

if ~has_Raman
    sp.Raman      = [];
    sp.Raman_cont = [];
end

% -------------------------------------------------------------------------
% Plot optional
% -------------------------------------------------------------------------
if do_plot
    local_plot_spectra(sp);
end

end  % G_spectra_nm


% =========================================================================
%  Local function for plotting (same rendering as G09_SPECTRA/G16_SPECTRA)
% =========================================================================
function local_plot_spectra(sp)

if ~isempty(sp.filename)
    [~, fname] = fileparts(sp.filename);
else
    fname = 'nm data';
end
fname_tex = strrep(fname, '_', '\_');

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

    for m = 1:sp.Nmodes
        if sp.Raman(m) > 0
            line(ax1, [sp.freq(m) sp.freq(m)], [0 sp.Raman(m)], ...
                 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, ...
                 'HandleVisibility', 'off');
        end
    end

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
