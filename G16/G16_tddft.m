function td = G16_tddft(filename, varargin)
% G16_TDDFT  Extracts TD-DFT excited states from a Gaussian 16 .out/.log file.
%
%   td = G16_TDDFT(filename)
%   td = G16_TDDFT(filename, 'nstates', N)     % load the first N states only
%   td = G16_TDDFT(filename, 'plot', true)     % also generate UV-Vis plot
%   td = G16_TDDFT(filename, 'FWHM_eV', 0.3)  % broadening FWHM (eV)
%
%   OUTPUT  struct td with fields:
%       .n           [Nstates×1]    excited state index (1, 2, ...)
%       .mult        {Nstates×1}    multiplicity label ('Singlet-A', 'Triplet-B',...)
%       .eV          [Nstates×1]    excitation energy (eV)
%       .nm          [Nstates×1]    excitation wavelength (nm)
%       .f           [Nstates×1]    oscillator strength (dimensionless)
%       .S2          [Nstates×1]    <S²> expectation value (0 for pure singlets)
%       .trans       {Nstates×1}    cell array: each entry is [Ncontrib×3]
%                                   columns: [MO_from, MO_to, coeff]
%                                   MO_to < 0 for de-excitation (←)
%       .Nstates     int            total number of excited states loaded
%       .has_S2      logical        true if <S²> values were found
%       .x_nm        [901×1]        wavelength grid in nm (100–1000 nm)
%       .eps_cont    [901×1]        Gaussian-broadened UV-Vis spectrum (arb. units)
%       .FWHM_eV     double         FWHM used for Gaussian broadening (eV)
%       .filename    char           source file path
%
%   Typical usage:
%       td = G16_tddft('violacein_td.out');
%       plot(td.x_nm, td.eps_cont)
%
%       % transition stem plot
%       stem(td.nm, td.f)
%
%       % transitions of the 2nd state
%       td.trans{2}   % [from, to, coeff]

% -------------------------------------------------------------------------
% Parse input arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename',  @ischar);
addParameter(p, 'nstates',   Inf,   @isnumeric);
addParameter(p, 'plot',      false, @islogical);
addParameter(p, 'FWHM_eV',  0.30,  @isnumeric);
parse(p, filename, varargin{:});

nstates_max = p.Results.nstates;
do_plot     = p.Results.plot;
FWHM_eV     = p.Results.FWHM_eV;

% -------------------------------------------------------------------------
% Read file
% -------------------------------------------------------------------------
if ~isfile(filename)
    error('G16_tddft: file not found: %s', filename);
end
fid  = fopen(filename,'r');
raw  = fread(fid,'*char')';
fclose(fid);
lines = strsplit(raw, newline);
G16_check_gaussian_match(lines, filename);
N = numel(lines);

% -------------------------------------------------------------------------
% Parsing stati eccitati
% -------------------------------------------------------------------------
% Formato:
%  Excited State   1:      Singlet-A      2.5661 eV  483.26 nm  f=0.0012  <S**2>=0.000
%       85 -> 86         0.69832
%       84 -> 86        -0.11234
%       85 <- 86         0.05123      (de-eccitazione, segno flip su to)

st_n    = [];
st_mult = {};
st_eV   = [];
st_nm   = [];
st_f    = [];
st_S2   = [];
st_trans = {};

k = 1;
while k <= N
    ln = lines{k};

    % ── "Excited State N: ..." header line
    % NOTE: <S**2>= used to be captured as a trailing optional group
    % (?:...)? tacked onto this same pattern, but MATLAB's regexp engine
    % silently drops a trailing optional capture group even when it
    % matches (confirmed with 'tokens','once', plain 'tokens', and
    % 'names' — all three drop it). Extracting it with its own
    % independent regexp below is the reliable fix.
    tok = regexp(ln, ...
        'Excited State\s+(\d+):\s+(\S+)\s+([\d.]+)\s+eV\s+([\d.]+)\s+nm\s+f=([\d.]+)', ...
        'tokens', 'once');

    if ~isempty(tok)
        sn   = str2double(tok{1});
        if sn > nstates_max
            k = k+1; continue
        end

        st_n(end+1)    = sn;                          %#ok<AGROW>
        st_mult{end+1} = tok{2};                      %#ok<AGROW>
        st_eV(end+1)   = str2double(tok{3});          %#ok<AGROW>
        st_nm(end+1)   = str2double(tok{4});          %#ok<AGROW>
        st_f(end+1)    = str2double(tok{5});          %#ok<AGROW>
        tokS2 = regexp(ln, '<S\*\*2>=([\d.]+)', 'tokens', 'once');
        if ~isempty(tokS2)
            st_S2(end+1) = str2double(tokS2{1});      %#ok<AGROW>
        else
            st_S2(end+1) = NaN;                       %#ok<AGROW>
        end

        % Read MO transition lines until a non-transition line is found
        trans_block = [];
        k = k+1;
        while k <= N
            ln2 = lines{k};
            % Transizione eccitante:  "  85 -> 86   0.69832"
            tk_ex = regexp(ln2, '^\s*(\d+)\s*->\s*(\d+)\s+([-\d.]+)', 'tokens', 'once');
            % De-eccitazione:         "  85 <- 86   0.05123"
            tk_de = regexp(ln2, '^\s*(\d+)\s*<-\s*(\d+)\s+([-\d.]+)', 'tokens', 'once');

            if ~isempty(tk_ex)
                from = str2double(tk_ex{1});
                to   = str2double(tk_ex{2});
                c    = str2double(tk_ex{3});
                trans_block(end+1,:) = [from, to, c]; %#ok<AGROW>
                k = k+1;
            elseif ~isempty(tk_de)
                from = str2double(tk_de{1});
                to   = -str2double(tk_de{2});   % negative = de-excitation
                c    = str2double(tk_de{3});
                trans_block(end+1,:) = [from, to, c]; %#ok<AGROW>
                k = k+1;
            else
                break   % end of transitions for this state
            end
        end
        st_trans{end+1} = trans_block;  %#ok<AGROW>
        continue
    end

    k = k+1;
end

% -------------------------------------------------------------------------
% Check
% -------------------------------------------------------------------------
Nstates = numel(st_n);
if Nstates == 0
    error('G16_tddft: no excited states found in %s\n(check that the job uses TD-DFT)', filename);
end

% -------------------------------------------------------------------------
% Spettro UV-Vis continuum (Gaussiana in energia, riportata su nm)
% -------------------------------------------------------------------------
% Wavelength grid from 100 to 1000 nm
x_nm   = (100 : 1 : 1000)';
x_eV_grid = 1239.84193 ./ x_nm;   % conversion nm -> eV (hc = 1239.84 eV·nm)

sigma_eV = FWHM_eV / (2*sqrt(2*log(2)));   % FWHM -> Gaussian sigma

eps_cont = zeros(size(x_nm));
for s = 1:Nstates
    if st_f(s) > 0
        eps_cont = eps_cont + st_f(s) * ...
            exp( -(x_eV_grid - st_eV(s)).^2 / (2*sigma_eV^2) );
    end
end

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
td.n         = st_n(:);
td.mult      = st_mult(:);
td.eV        = st_eV(:);
td.nm        = st_nm(:);
td.f         = st_f(:);
td.S2        = st_S2(:);
td.trans     = st_trans(:);
td.Nstates   = Nstates;
td.has_S2    = ~all(isnan(st_S2));
td.x_nm      = x_nm;
td.eps_cont  = eps_cont;
td.FWHM_eV   = FWHM_eV;
td.filename  = filename;

% -------------------------------------------------------------------------
% Print summary
% -------------------------------------------------------------------------
fprintf('\n── G16_tddft: %s ──\n', filename);
fprintf('  %d excited states read\n', Nstates);
fprintf('  %-6s  %-14s  %8s  %8s  %8s\n', 'State', 'Mult', 'eV', 'nm', 'f');
fprintf('  %s\n', repmat('-',1,52));
for s = 1:min(Nstates, 20)
    fprintf('  %-6d  %-14s  %8.4f  %8.2f  %8.4f\n', ...
        td.n(s), td.mult{s}, td.eV(s), td.nm(s), td.f(s));
end
if Nstates > 20
    fprintf('  ... (%d states total)\n', Nstates);
end
fprintf('\n');

% -------------------------------------------------------------------------
% Plot optional
% -------------------------------------------------------------------------
if do_plot
    [~, fname] = fileparts(filename);

    fig = figure('Color','white','Name',fname,'NumberTitle','off');
    ax  = axes('Parent', fig);
    hold(ax, 'on');

    % stick f
    for s = 1:Nstates
        if td.f(s) > 0
            line(ax, [td.nm(s) td.nm(s)], [0, td.f(s)], ...
                 'Color', [0.70 0.70 0.70], 'LineWidth', 0.9, ...
                 'HandleVisibility', 'off');
        end
    end

    % continuum scaled to f_max
    f_max   = max(td.f);
    eps_max = max(td.eps_cont);
    if eps_max > 0
        scale = f_max / eps_max;
    else
        scale = 1;
    end
    plot(ax, td.x_nm, td.eps_cont * scale, ...
         'Color', [0.15 0.45 0.80], 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Gaussiana FWHM=%.2f eV', FWHM_eV));

    set(ax, 'Box','on', 'XGrid','off', 'YGrid','off');
    xlabel(ax, 'Wavelength (nm)', 'FontSize',10);
    ylabel(ax, 'Oscillator strength  f', 'FontSize',10);
    title(ax,  strrep(fname,'_','\_'), 'FontSize',11, 'Interpreter','tex');
    legend(ax, 'show', 'Location','northeast', 'Box','off');
    xlim(ax, [200, max(td.nm)*1.15]);
end

end  % G16_tddft
