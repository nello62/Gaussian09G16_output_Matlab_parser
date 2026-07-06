function data = G09_fchk_read(filename, varargin)
% G09_FCHK_READ  Reads a Gaussian 09/16 formatted checkpoint file (.fchk).
%
%   data = G09_FCHK_READ(filename)
%   data = G09_FCHK_READ(filename, 'verbose', false)
%
%   The .fchk file is the ASCII-formatted version of the binary .chk file,
%   generated with the Gaussian utility:
%       formchk  jobname.chk  jobname.fchk
%
%   DIFFERENCE FROM g09_fckread.m (old .fck reader):
%     - Reads the standard .fchk format (not the obsolete .fck / newzmat format)
%     - Uses a single generic parser for all sections (no hard-coded line counts)
%     - Returns SI-ready fields: XYZ in Angstrom, dipole in Debye
%     - No external helper functions required (self-contained)
%     - Compatible with both G09 and G16 .fchk files
%
%   Optional parameters (Name-Value):
%       'verbose'  - true (default) | false  — print progress messages
%
%   OUTPUT  struct data with fields:
%
%   ── Molecular information ─────────────────────────────────────────────
%       .title        char          first line of the .fchk file
%       .method       char          method string (e.g. 'RB3LYP')
%       .basis        char          basis set string (e.g. 'def2SVPP')
%       .Nat          int           number of atoms
%       .charge       int           molecular charge
%       .mult         int           spin multiplicity
%       .Nelec        int           total number of electrons
%       .Nalpha       int           number of alpha electrons
%       .Nbeta        int           number of beta electrons
%       .Nbasis       int           number of basis functions
%       .Nbasis_indep int           number of independent basis functions
%
%   ── Geometry ──────────────────────────────────────────────────────────
%       .symbols      {Nat x 1}     atomic symbols (e.g. 'C', 'H', 'N')
%       .AN           [Nat x 1]     atomic numbers
%       .masses       [Nat x 1]     real atomic masses (AMU)
%       .xyz          [Nat x 3]     Cartesian coordinates (Angstrom)
%       .xyz_bohr     [Nat x 3]     Cartesian coordinates (Bohr)
%
%   ── Energies ──────────────────────────────────────────────────────────
%       .SCF_energy   double        SCF energy (Hartree)
%       .total_energy double        total energy (Hartree)
%       .virial_ratio double        virial ratio (should be ~2)
%
%   ── Electronic structure ──────────────────────────────────────────────
%       .alpha_orb_energies  [Nbasis x 1]   alpha MO energies (Hartree)
%       .beta_orb_energies   [Nbasis x 1]   beta MO energies  (Hartree), [] if RHF
%       .alpha_MO_coeff      [Nbasis^2 x 1] alpha MO coefficients (column-major)
%       .mulliken_charges    [Nat x 1]      Mulliken charges (e)
%       .HOMO_idx            int            index of HOMO (1-based)
%       .HOMO_eV             double         HOMO energy (eV)
%       .LUMO_eV             double         LUMO energy (eV)
%       .gap_eV              double         HOMO-LUMO gap (eV)
%
%   ── Forces and geometry optimisation ─────────────────────────────────
%       .gradient     [Nat x 3]     Cartesian gradient (Hartree/Bohr)
%       .rms_force    double        RMS force
%       .force_const  [3Nat x 3Nat] Cartesian force constant matrix (Hartree/Bohr^2)
%                                   reconstructed from lower triangle
%
%   ── Properties ────────────────────────────────────────────────────────
%       .dipole_au    [1 x 3]       dipole moment (atomic units, e·a0)
%       .dipole_D     [1 x 3]       dipole moment (Debye)
%       .dipole_tot_D double        dipole magnitude (Debye)
%
%       .polar_au     [3 x 3]       polarisability tensor α (au, symmetric)
%       .polar_iso    double        isotropic polarisability (au)
%       .polar_aniso  double        polarisability anisotropy (au)
%
%       .beta_au      struct        first hyperpolarisability β (au)
%                       .xxx .xxy .xyy .yyy .xxz .xyz .yyz .xzz .yzz .zzz
%       .beta_vec     double        |β_vec| = sqrt(βx²+βy²+βz²) (au)
%
%       .dipole_deriv [3 x 3Nat]    dipole derivatives dμ/dR (au, [x/y/z, 3Nat])
%                                   dmu_x = data.dipole_deriv(1,:)
%
%       .polar_deriv  [6 x 3Nat]    polarisability derivatives dα/dR (au)
%                                   rows: xx,xy,yy,xz,yz,zz
%
%   .filename   char                source file path
%
%   Example:
%       data = G09_fchk_read('3typ.fchk');
%       data.gap_eV                     % HOMO-LUMO gap in eV
%       data.dipole_D                   % dipole vector in Debye
%       data.polar_iso                  % isotropic polarisability in au
%       data.force_const(1:6,1:6)       % top-left of force constant matrix
%       data.mulliken_charges           % Mulliken charges

% -------------------------------------------------------------------------
% Parse arguments
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p,  'filename', @ischar);
addParameter(p, 'verbose',  true,  @islogical);
parse(p, filename, varargin{:});
verbose = p.Results.verbose;

if ~isfile(filename)
    error('G09_fchk_read: file not found: %s', filename);
end

% -------------------------------------------------------------------------
% Z -> symbol table
% -------------------------------------------------------------------------
SYM = { ...
    'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', ...
    'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca', ...
    'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', ...
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr', ...
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', ...
    'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', ...
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', ...
    'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', ...
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', ...
    'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', ...
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', ...
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'};

% -------------------------------------------------------------------------
% Read file
% -------------------------------------------------------------------------
fid   = fopen(filename, 'r');
raw   = fread(fid, '*char')';
fclose(fid);
lines = strsplit(raw, newline);
N     = numel(lines);

% -------------------------------------------------------------------------
% Generic .fchk section parser
% -------------------------------------------------------------------------
% .fchk format:
%   Scalar:  "<label padded to 40 chars>  I|R|C  <value>"
%   Array:   "<label padded to 40 chars>  I|R|C  N=  <count>"
%            followed by data lines (5 values/line for R, 6 for I)

% Build index of all section header lines
%   sec_map(name) = {line_index, type, n_values}
sec_lines   = zeros(1, N, 'logical');
sec_names   = {};
sec_types   = {};    % 'I', 'R', 'C'
sec_counts  = [];    % 0 for scalar
sec_idx     = [];    % line index in 'lines' (1-based)

header_re = '^(.{1,43}?)\s{1,3}(I|R|C)\s+(N=\s*(\d+)|[-\d.E+]+)';

for k = 1:N
    ln = lines{k};
    if numel(ln) < 45, continue; end
    tok = regexp(ln, header_re, 'tokens', 'once');
    if isempty(tok), continue; end
    name  = strtrim(tok{1});
    dtype = tok{2};
    rest  = strtrim(tok{3});
    if startsWith(rest, 'N=')
        nvals = str2double(strtrim(rest(3:end)));
    else
        nvals = 0;   % scalar — value is on the same line
    end
    sec_names{end+1}  = name;  %#ok<AGROW>
    sec_types{end+1}  = dtype; %#ok<AGROW>
    sec_counts(end+1) = nvals; %#ok<AGROW>
    sec_idx(end+1)    = k;     %#ok<AGROW>
end

% Helper: read a named section
    function val = read_sec(keyword)
        % Find the LAST matching header (handles duplicate sections)
        idx_match = find(~cellfun(@isempty, ...
            regexpi(sec_names, ['^', regexptranslate('wildcard', keyword), '$'])));
        if isempty(idx_match)
            % Fallback: substring match
            idx_match = find(~cellfun(@isempty, ...
                cellfun(@(s) strfind(lower(s), lower(keyword)), sec_names, 'UniformOutput', false)));
        end
        if isempty(idx_match)
            val = [];
            return
        end
        mi    = idx_match(end);
        k0    = sec_idx(mi);
        dtype = sec_types{mi};
        nvals = sec_counts(mi);

        if nvals == 0
            % Scalar: parse from the header line itself
            tok = regexp(lines{k0}, ...
                '(I|R|C)\s+([-\d.E+]+)\s*$', 'tokens', 'once');
            if isempty(tok)
                val = [];
            elseif strcmp(dtype, 'I')
                val = str2double(tok{2});
            else
                val = str2double(tok{2});
            end
            return
        end

        % Array: read data lines
        vals = [];
        k2   = k0 + 1;
        while numel(vals) < nvals && k2 <= N
            ln2 = strtrim(lines{k2});
            if isempty(ln2), k2 = k2+1; continue; end
            % Stop if next header line reached
            if numel(ln2) > 44 && ~isempty(regexp(ln2, header_re, 'once'))
                break
            end
            v = sscanf(ln2, '%f');
            vals = [vals; v]; %#ok<AGROW>
            k2 = k2 + 1;
        end
        val = vals(1:min(end, nvals));
    end

% -------------------------------------------------------------------------
% Extract header info (lines 1-2)
% -------------------------------------------------------------------------
title_line = strtrim(lines{1});
hdr2       = strtrim(lines{2});
% Line 2: "CalcType  Method  BasisSet"
tok2 = strsplit(hdr2);
if numel(tok2) >= 3
    calc_type = tok2{1};
    method    = tok2{2};
    basis     = strjoin(tok2(3:end), ' ');
elseif numel(tok2) == 2
    method = tok2{1}; basis = tok2{2}; calc_type = '';
else
    method = hdr2; basis = ''; calc_type = '';
end

% -------------------------------------------------------------------------
% Scalar quantities
% -------------------------------------------------------------------------
Nat          = round(read_sec('Number of atoms'));
charge       = round(read_sec('Charge'));
mult         = round(read_sec('Multiplicity'));
Nelec        = round(read_sec('Number of electrons'));
Nalpha       = round(read_sec('Number of alpha electrons'));
Nbeta        = round(read_sec('Number of beta electrons'));
Nbasis       = round(read_sec('Number of basis functions'));
Nbasis_indep = round(read_sec('Number of independent functions'));
SCF_energy   = read_sec('SCF Energy');
total_energy = read_sec('Total Energy');
virial       = read_sec('Virial Ratio');
rms_force    = read_sec('RMS Force');

if verbose
    fprintf('\n── G09_fchk_read: %s ──\n', filename);
    fprintf('  Title  : %s\n', title_line);
    fprintf('  Method : %s  Basis: %s\n', method, basis);
    fprintf('  Nat=%d  Charge=%+d  Mult=%d  Nbasis=%d\n', Nat, charge, mult, Nbasis);
    fprintf('  SCF Energy = %.10f Ha\n', SCF_energy);
end

% -------------------------------------------------------------------------
% Geometry
% -------------------------------------------------------------------------
AN_vec   = round(read_sec('Atomic numbers'));
masses   = read_sec('Real atomic weights');
xyz_bohr = read_sec('Current cartesian coordinates');

if numel(AN_vec) ~= Nat
    warning('G09_fchk_read: Atomic numbers count mismatch.');
end

% Convert Bohr -> Angstrom (CODATA 2022 Bohr radius)
a0_ang   = 0.529177210544;
xyz_bohr = reshape(xyz_bohr, 3, Nat)';    % (Nat x 3)
xyz_ang  = xyz_bohr * a0_ang;

% Build symbol list
symbols = cell(Nat, 1);
for i = 1:Nat
    z = AN_vec(i);
    if z >= 1 && z <= numel(SYM)
        symbols{i} = SYM{z};
    else
        symbols{i} = sprintf('Z%d', z);
    end
end

if verbose
    fprintf('  Atoms  :');
    for i = 1:min(Nat, 8)
        fprintf(' %s', symbols{i});
    end
    if Nat > 8, fprintf(' ... (%d total)', Nat); end
    fprintf('\n');
end

% -------------------------------------------------------------------------
% Electronic structure
% -------------------------------------------------------------------------
alpha_orb = read_sec('Alpha Orbital Energies');
beta_orb  = read_sec('Beta Orbital Energies');    % [] for RHF
alpha_MO  = read_sec('Alpha MO coefficients');

mull_chg  = read_sec('Mulliken Charges');

% HOMO / LUMO
HOMO_idx = Nalpha;   % 1-based index of HOMO
HOMO_eV  = NaN; LUMO_eV = NaN; gap_eV = NaN;
if numel(alpha_orb) >= HOMO_idx + 1
    Eh2eV   = 27.211386;
    HOMO_eV = alpha_orb(HOMO_idx)     * Eh2eV;
    LUMO_eV = alpha_orb(HOMO_idx + 1) * Eh2eV;
    gap_eV  = LUMO_eV - HOMO_eV;
    if verbose
        fprintf('  HOMO   = %.4f eV   LUMO = %.4f eV   Gap = %.4f eV\n', ...
            HOMO_eV, LUMO_eV, gap_eV);
    end
end

% -------------------------------------------------------------------------
% Forces
% -------------------------------------------------------------------------
grad_raw  = read_sec('Cartesian Gradient');
gradient  = NaN(Nat, 3);
if numel(grad_raw) == 3*Nat
    gradient = reshape(grad_raw, 3, Nat)';   % (Nat x 3) in Hartree/Bohr
end

% Force constant matrix: lower triangle, ROW by row (Gaussian .fchk convention)
% Storage order: F(1,1), F(2,1), F(2,2), F(3,1), F(3,2), F(3,3), ...
% i.e. for row i (1..N3): elements F(i,1), F(i,2), ..., F(i,i)
fc_raw = read_sec('Cartesian Force Constants');
N3     = 3 * Nat;
force_const = NaN(N3, N3);
if numel(fc_raw) == N3*(N3+1)/2
    k_fc = 1;
    for row = 1:N3
        for col = 1:row          % lower triangle: col ≤ row
            force_const(row, col) = fc_raw(k_fc);
            force_const(col, row) = fc_raw(k_fc);   % symmetric
            k_fc = k_fc + 1;
        end
    end
end

% -------------------------------------------------------------------------
% Dipole moment
% -------------------------------------------------------------------------
dip_raw   = read_sec('Dipole Moment');       % [x y z] in au (e·a0)
dip_au    = [NaN NaN NaN];
dip_D     = [NaN NaN NaN];
dip_tot_D = NaN;

if numel(dip_raw) == 3
    au2D     = 2.541747;
    dip_au   = dip_raw(:)';
    dip_D    = dip_au * au2D;
    dip_tot_D = norm(dip_D);
    if verbose
        fprintf('  Dipole  = (%.4f, %.4f, %.4f) D   |μ| = %.4f D\n', ...
            dip_D(1), dip_D(2), dip_D(3), dip_tot_D);
    end
end

% -------------------------------------------------------------------------
% Polarisability
% -------------------------------------------------------------------------
% .fchk stores upper triangle by rows: xx xy yy xz yz zz
pol_raw  = read_sec('Polarizability');
pol_au   = NaN(3,3);
pol_iso  = NaN;
pol_aniso = NaN;

if numel(pol_raw) == 6
    v = pol_raw;
    % Order: xx(1) xy(2) yy(3) xz(4) yz(5) zz(6)
    pol_au = [v(1) v(2) v(4); ...
              v(2) v(3) v(5); ...
              v(4) v(5) v(6)];
    pol_iso   = trace(pol_au) / 3;
    pol_aniso = sqrt(0.5 * ( ...
        (pol_au(1,1)-pol_au(2,2))^2 + ...
        (pol_au(2,2)-pol_au(3,3))^2 + ...
        (pol_au(3,3)-pol_au(1,1))^2 + ...
        6*(pol_au(1,2)^2 + pol_au(1,3)^2 + pol_au(2,3)^2)));
    if verbose
        fprintf('  α_iso   = %.4f au   α_aniso = %.4f au\n', pol_iso, pol_aniso);
    end
end

% -------------------------------------------------------------------------
% First hyperpolarisability
% -------------------------------------------------------------------------
% .fchk stores 10 unique components: xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz
hyper_raw = read_sec('HyperPolarizability');
beta      = struct('xxx',NaN,'xxy',NaN,'xyy',NaN,'yyy',NaN,'xxz',NaN, ...
                   'xyz',NaN,'yyz',NaN,'xzz',NaN,'yzz',NaN,'zzz',NaN);
beta_vec  = NaN;

if numel(hyper_raw) >= 10
    beta.xxx = hyper_raw(1);
    beta.xxy = hyper_raw(2);
    beta.xyy = hyper_raw(3);
    beta.yyy = hyper_raw(4);
    beta.xxz = hyper_raw(5);
    beta.xyz = hyper_raw(6);
    beta.yyz = hyper_raw(7);
    beta.xzz = hyper_raw(8);
    beta.yzz = hyper_raw(9);
    beta.zzz = hyper_raw(10);
    % Vector component: β_i = (1/5) Σ_j (β_ijj + β_jij + β_jji)
    % For parallel component β_x,y,z:
    bx = (beta.xxx + beta.xyy + beta.xzz);
    by = (beta.xxy + beta.yyy + beta.yzz);
    bz = (beta.xxz + beta.yyz + beta.zzz);
    beta_vec = sqrt(bx^2 + by^2 + bz^2);
    if verbose
        fprintf('  |β_vec| = %.2f au\n', beta_vec);
    end
end

% -------------------------------------------------------------------------
% Dipole derivatives  dμ/dR  (APT tensor)
% -------------------------------------------------------------------------
% Stored as 3x(3Nat) values: [dμx/dR1x, dμy/dR1x, dμz/dR1x, dμx/dR1y, ...]
% i.e. for each Cartesian displacement, the 3 dipole derivative components
dip_deriv_raw = read_sec('Dipole Derivatives');
dip_deriv     = NaN(3, N3);   % rows: x,y,z  cols: 3Nat displacements

if numel(dip_deriv_raw) == 3*N3
    % Column order in fchk: for each displacement k, 3 components (μx,μy,μz)
    tmp = reshape(dip_deriv_raw, 3, N3);   % (3 x 3Nat)
    dip_deriv = tmp;                        % row i = dμ_i/dR_j for all j
end

% -------------------------------------------------------------------------
% Polarisability derivatives  dα/dR
% -------------------------------------------------------------------------
% Stored as 6x(3Nat) values: 6 unique α components per displacement
% Order per displacement: xx,xy,yy,xz,yz,zz
pol_deriv_raw = read_sec('Polarizability Derivatives');
pol_deriv     = NaN(6, N3);   % rows: xx,xy,yy,xz,yz,zz

if numel(pol_deriv_raw) == 6*N3
    pol_deriv = reshape(pol_deriv_raw, 6, N3);
end

% -------------------------------------------------------------------------
% Assemble output struct
% -------------------------------------------------------------------------
data.title        = title_line;
data.method       = method;
data.basis        = basis;
data.calc_type    = calc_type;
data.Nat          = Nat;
data.charge       = charge;
data.mult         = mult;
data.Nelec        = Nelec;
data.Nalpha       = Nalpha;
data.Nbeta        = Nbeta;
data.Nbasis       = Nbasis;
data.Nbasis_indep = Nbasis_indep;

data.symbols      = symbols;
data.AN           = AN_vec;
data.masses       = masses;
data.xyz          = xyz_ang;
data.xyz_bohr     = xyz_bohr;

data.SCF_energy   = SCF_energy;
data.total_energy = total_energy;
data.virial_ratio = virial;
data.rms_force    = rms_force;

data.alpha_orb_energies = alpha_orb;
data.beta_orb_energies  = beta_orb;
data.alpha_MO_coeff     = alpha_MO;
data.mulliken_charges   = mull_chg;
data.HOMO_idx           = HOMO_idx;
data.HOMO_eV            = HOMO_eV;
data.LUMO_eV            = LUMO_eV;
data.gap_eV             = gap_eV;

data.gradient     = gradient;
data.force_const  = force_const;

data.dipole_au    = dip_au;
data.dipole_D     = dip_D;
data.dipole_tot_D = dip_tot_D;

data.polar_au     = pol_au;
data.polar_iso    = pol_iso;
data.polar_aniso  = pol_aniso;

data.beta_au      = beta;
data.beta_vec     = beta_vec;

data.dipole_deriv = dip_deriv;
data.polar_deriv  = pol_deriv;

data.filename     = filename;

% -------------------------------------------------------------------------
% Build compatibility sub-structs for G09_draw_molecule, G09_charges,
% G09_draw_mode
% -------------------------------------------------------------------------

% ── mol  (compatible with G09_draw_molecule, G09_draw_mode) ──────────────
mol.symbols     = symbols;
mol.xyz         = xyz_ang;
mol.Z           = AN_vec;
mol.Natoms      = Nat;
mol.step        = 1;
mol.n_steps     = 1;
mol.orientation = 'fchk (Input orientation)';
mol.filename    = filename;
data.mol        = mol;

% ── ch  (compatible with G09_charges) ────────────────────────────────────
% Mulliken charges are always in the .fchk file
ch.symbols   = symbols;
ch.charges   = mull_chg;
ch.charges_H = [];          % H-summed not in fchk — use G09_charges on .log instead
ch.sum_q     = sum(mull_chg);
ch.type      = 'Mulliken';
ch.label     = 'Mulliken Charges (from .fchk)';
ch.Natoms    = Nat;
ch.filename  = filename;
data.ch      = ch;

% ── nm  (compatible with G09_draw_mode) ───────────────────────────────────
% Compute normal modes from the mass-weighted Cartesian force constant matrix.
% The .fchk force_const is in Hartree/Bohr^2.
%
% Steps:
%   1. Build mass-weight matrix M^{-1/2} (diagonal, 3Nat x 3Nat)
%   2. Mass-weight the Hessian: Fmw = M^{-1/2} * F * M^{-1/2}
%   3. Diagonalise Fmw: eigenvalues λ_i (Hartree/(Bohr^2·amu_e units)
%   4. Frequencies: ω_i = sqrt(λ_i) in au  →  cm⁻¹
%   5. Eigenvectors: displacement vectors (mass-weighted)  →  Cartesian

nm = struct();
nm.Nmodes    = 0;
nm.Natoms    = Nat;
nm.has_Raman = false;
nm.freq      = [];
nm.IR        = [];
nm.Raman     = [];
nm.disp      = [];
nm.filename  = filename;

if ~any(isnan(force_const(:))) && numel(masses) == Nat
    % ── Mass vector (amu) ────────────────────────────────────────────────
    m_vec = repelem(masses(:), 3);      % [3Nat x 1] in amu

    % ── Build translation/rotation (TR) projector ────────────────────────
    % Wilson-Decius-Cross procedure:
    %   D columns = mass-weighted translation and rotation vectors
    %   P_vib = I - D_orth * D_orth'  projects onto the vibrational subspace
    % Coordinates must be in Bohr, centred at centre of mass.
    M_tot  = sum(masses);
    com    = (masses(:)' * xyz_bohr) / M_tot;     % [1x3] Bohr
    xyz_c  = xyz_bohr - com;                       % centred coordinates

    D_tr = zeros(N3, 6);
    for ii = 1:Nat
        sqm = sqrt(masses(ii));
        xi = xyz_c(ii,1); yi = xyz_c(ii,2); zi = xyz_c(ii,3);
        D_tr(3*ii-2, 1) = sqm;              % Tx
        D_tr(3*ii-1, 2) = sqm;              % Ty
        D_tr(3*ii,   3) = sqm;              % Tz
        D_tr(3*ii-2, 4) = 0;               % Rx: (0, -z,  y) × sqrt(m)
        D_tr(3*ii-1, 4) = -zi * sqm;
        D_tr(3*ii,   4) =  yi * sqm;
        D_tr(3*ii-2, 5) =  zi * sqm;       % Ry: ( z,  0, -x) × sqrt(m)
        D_tr(3*ii-1, 5) = 0;
        D_tr(3*ii,   5) = -xi * sqm;
        D_tr(3*ii-2, 6) = -yi * sqm;       % Rz: (-y,  x,  0) × sqrt(m)
        D_tr(3*ii-1, 6) =  xi * sqm;
        D_tr(3*ii,   6) = 0;
    end

    % Orthonormalise (QR) — removes linearly dependent columns,
    % e.g. for linear molecules only 2 rotation vectors are non-zero
    [Q_tr, R_tr] = qr(D_tr, 0);
    tr_keep = abs(diag(R_tr)) > 1e-6;
    D_orth  = Q_tr(:, tr_keep);            % [3Nat × n_tr]
    n_tr    = size(D_orth, 2);             % 5 (linear) or 6 (non-linear)

    % Projector onto vibrational subspace
    P_vib = eye(N3) - D_orth * D_orth';

    % ── Mass-weighted, projected Hessian ─────────────────────────────────
    mw    = 1.0 ./ sqrt(m_vec);
    Fmw   = force_const .* (mw * mw');
    Fmw   = (Fmw + Fmw') / 2;

    Fmw_proj = P_vib * Fmw * P_vib;
    Fmw_proj = (Fmw_proj + Fmw_proj') / 2;

    % ── Diagonalise ───────────────────────────────────────────────────────
    [V, D_eig] = eig(Fmw_proj);
    lambda     = diag(D_eig);              % eigenvalues in Eh/(Bohr²·amu)

    % ── Convert to cm⁻¹ ─────────────────────────────────────────────────
    % ν (cm⁻¹) = sqrt(λ · Eh_J / (a0_m² · amu_kg)) / (2π · c_cms)
    % hess2cm1 ≈ 5140.49 cm⁻¹ per sqrt(Eh/(Bohr²·amu))
    Eh_J     = 4.3597447222060e-18;  % J        CODATA 2022
    a0_m     = 5.29177210544e-11;    % m
    amu_kg   = 1.66053906892e-27;    % kg
    c_cms    = 2.99792458e10;        % cm/s
    hess2cm1 = sqrt(Eh_J / (a0_m^2 * amu_kg)) / (2*pi * c_cms);

    freq_cm1 = sign(lambda) .* sqrt(abs(lambda)) * hess2cm1;

    % ── Sign convention note ──────────────────────────────────────────────
    % Eigenvector sign from eig() is numerically arbitrary (LAPACK does not
    % guarantee a fixed convention). A mode and its sign-flipped counterpart
    % represent the SAME physical vibration (180° phase difference), but
    % may appear reversed compared to Gaussian's own .out file or GaussView.
    % Use G09_draw_mode(..., 'FlipSign', true) to invert arrows for a
    % specific mode if it appears antiparallel to GaussView's rendering.

    % ── Cartesian displacement vectors ───────────────────────────────────
    % L_cart = M^{-1/2} * V  then normalised to Cartesian unit length,
    % matching Gaussian's print convention (sum L_i² = 1 per mode).
    % The .out file prints these vectors rounded to 2 decimal places;
    % agreement with .fchk-computed vectors is within ±0.005 (rounding only).
    L_cart    = mw .* V;
    col_norms = sqrt(sum(L_cart.^2, 1));
    col_norms(col_norms == 0) = 1;
    L_cart    = L_cart ./ col_norms;   % unit Cartesian norm — same as Gaussian

    % ── Sort and discard TR modes ─────────────────────────────────────────
    [freq_sorted, sort_idx] = sort(freq_cm1);
    L_sorted    = L_cart(:, sort_idx);

    % Skip the n_tr near-zero eigenvalue modes (translations + rotations)
    freq_vib    = freq_sorted(n_tr+1:end);
    L_vib_cn    = L_sorted(:, n_tr+1:end);  % Cartesian unit-norm (matches .out)
    Nmodes_vib  = N3 - n_tr;

    % Mass-normalised L for intensity calculations (sum m*L²=1 per mode)
    % Recompute from V before renormalisation: L_mn = mw .* V
    L_mn_sorted = mw .* V(:, sort_idx);     % (N3 × N3), NOT unit-normalised
    L_vib_mn    = L_mn_sorted(:, n_tr+1:end);

    % Reshape displacements: [Nat × 3 × Nmodes_vib]  — unit-norm convention
    disp_all    = reshape(L_vib_cn, 3, Nat, Nmodes_vib);
    disp_all    = permute(disp_all, [2, 1, 3]);

    % ── IR intensities (KM/Mole) ─────────────────────────────────────────
    % Storage order in .fchk: (dμx/dRk, dμy/dRk, dμz/dRk) for k=1..3Nat
    %   → reshape(N3,3) then transpose → (3 × N3).
    % L_vib_mn is mass-normalised (sum m*L²=1) → dμ/dQ in e/sqrt(amu).
    % Conversion: C_IR = 42.2561 * (2.541747/0.529177)² = 974.88 KM/mol
    %             per (e/sqrt(amu))²
    C_IR   = 42.2561 * (2.541747 / 0.529177)^2;
    IR_int = zeros(Nmodes_vib, 1);
    if ~any(isnan(dip_deriv(:)))
        DD     = reshape(dip_deriv, N3, 3)';   % (3 × N3): correct storage
        dmu_dQ = DD * L_vib_mn;                % (3 × Nmodes_vib)
        IR_int = sum(dmu_dQ.^2, 1)' * C_IR;
    end

    % ── Raman activities (Å⁴/AMU) ────────────────────────────────────────
    % Storage order in .fchk: (αxx,αxy,αyy,αxz,αyz,αzz) for k=1..3Nat
    %   → reshape(N3,6) then transpose → (6 × N3).
    % Conversion: C_Ra = (0.148185/0.529177)² = 0.07841 Å⁴/amu per au²/amu
    C_Ra      = (0.148185 / 0.529177)^2;
    Raman_int = zeros(Nmodes_vib, 1);
    has_Raman = false;
    if ~any(isnan(pol_deriv(:)))
        PD     = reshape(pol_deriv, N3, 6)';   % (6 × N3): correct storage
        da_dQ  = PD * L_vib_mn;                % (6 × Nmodes_vib)
        d_iso  = (da_dQ(1,:) + da_dQ(3,:) + da_dQ(6,:)) / 3;
        d_an2  = ((da_dQ(1,:)-da_dQ(3,:)).^2 + ...
                  (da_dQ(3,:)-da_dQ(6,:)).^2 + ...
                  (da_dQ(6,:)-da_dQ(1,:)).^2) / 2 + ...
                 3*(da_dQ(2,:).^2 + da_dQ(4,:).^2 + da_dQ(5,:).^2);
        Raman_int = (45*d_iso.^2 + 7*d_an2)' * C_Ra;
        has_Raman = true;
    end

    % Fill nm struct — same fields as G09_nmodes output
    nm.Nmodes    = Nmodes_vib;
    nm.Natoms    = Nat;
    nm.has_Raman = has_Raman;
    nm.freq      = freq_vib;
    nm.IR        = IR_int;
    nm.Raman     = Raman_int;
    nm.disp      = disp_all;    % [Nat x 3 x Nmodes_vib]
    nm.symmetry  = repmat({'?'}, Nmodes_vib, 1);
    nm.redmass   = [];
    nm.frcconst  = [];
    nm.filename  = filename;

    if verbose
        fprintf('  Normal modes: %d vibrational  (%d TR projected out)\n', ...
            Nmodes_vib, n_tr);
        first_real = find(freq_vib > 10, 1);
        if ~isempty(first_real)
            fprintf('  Lowest freq  : %.1f cm-1 (mode %d)\n', ...
                freq_vib(first_real), first_real);
            fprintf('  Highest freq : %.1f cm-1\n', freq_vib(end));
        end
        n_imag = sum(freq_vib < -5);
        if n_imag > 0
            fprintf('  Imaginary (<-5 cm-1): %d  (saddle point / TS geometry)\n', n_imag);
        end
    end
end

data.nm = nm;

if verbose
    fprintf('\n  Compatibility sub-structs ready:\n');
    fprintf('    data.mol  → G09_draw_molecule(data.mol)\n');
    fprintf('    data.ch   → G09_charges uses data.ch directly\n');
    fprintf('    data.nm   → G09_draw_mode(data.mol, data.nm, mode_idx)\n');
    fprintf('  Done.\n\n');
end

end  % G09_fchk_read
