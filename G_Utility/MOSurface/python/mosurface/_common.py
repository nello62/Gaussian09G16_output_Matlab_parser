"""Shared low-level helpers for the MOSurface Python port: the Boys
function, a percentile helper, isosurface-clipping/shell-support checks,
and Gaussian .cube/.fchk-basis-set file I/O. Not part of the public API
-- import the four g16_draw_*_surface functions from the top-level
``mosurface`` package instead.

Self-contained by design, mirroring the MATLAB G_Utility/MOSurface
philosophy: no dependency on G16parser or any core-toolbox package, so
this companion can be copied/used independently of either.
"""

import re

import numpy as np
from scipy.special import gammainc, gammaln


def boys_function(x, nmax):
    """Boys function F_n(x) = int_0^1 t^(2n) exp(-x t^2) dt for n=0..nmax,
    via the exact closed form in terms of the regularized lower
    incomplete gamma function:
        F_n(x) = gammainc(n+1/2, x) * Gamma(n+1/2) / (2 x^(n+1/2)),  x > 0
        F_n(0) = 1/(2n+1)
    Note scipy.special.gammainc(a, x) takes the shape parameter FIRST,
    unlike MATLAB's gammainc(x, a) -- this is the most common porting
    slip for this function, so it is called out explicitly here.

    x : array-like, any shape. Returns F: [x.size, nmax+1], column n =
    F_n(x), reshaped back to x's own shape with an extra trailing axis.
    """
    x = np.asarray(x, dtype=float).ravel()
    npts = x.size
    F = np.zeros((npts, nmax + 1))
    small = x < 1e-12

    for n in range(nmax + 1):
        col = np.zeros(npts)
        a = n + 0.5
        if np.any(~small):
            xv = x[~small]
            # Gamma(a) = exp(gammaln(a)); avoids scipy.special.gamma
            # overflow for the (here unneeded) large-a case, and matches
            # MATLAB's gammainc(x,a)*gamma(a)/(2*x^a) construction exactly.
            col[~small] = gammainc(a, xv) * np.exp(gammaln(a)) / (2 * xv**a)
        if np.any(small):
            col[small] = 1.0 / (2 * n + 1)
        F[:, n] = col
    return F


def pctile_local(v, p):
    """P-th percentile of v via a plain sort (same nearest-rank indexing
    convention as MATLAB's g_pctile_local, NOT numpy's default linear
    interpolation -- kept for numerical consistency with the MATLAB
    original rather than for any correctness reason).
    """
    v = np.sort(np.asarray(v).ravel())
    n = v.size
    if n == 0:
        return 0.0
    idx = int(np.clip(round(p / 100.0 * n), 1, n)) - 1
    return v[idx]


def is_clipped(vertices, Xa, Ya, Za, tol):
    """True if any isosurface vertex lies within tol of the evaluated
    grid's outer boundary -- a sign the true isosurface extends beyond
    the grid and has been artificially cut off there.
    """
    if vertices is None or len(vertices) == 0:
        return False
    lo = np.array([Xa.min(), Ya.min(), Za.min()])
    hi = np.array([Xa.max(), Ya.max(), Za.max()])
    v = np.asarray(vertices)
    return bool(np.any((v <= lo + tol) | (v >= hi - tol)))


def check_supported_shells(shell_types):
    """Raises a clear error (no partial rendering) if the basis contains
    any shell type beyond S/P/SP/D/F/G (pure or Cartesian).
    """
    supported = {0, 1, -1, -2, 2, -3, 3, -4, 4}
    bad = sorted({int(t) for t in shell_types if int(t) not in supported})
    if not bad:
        return
    raise ValueError(
        f"check_supported_shells: this basis set contains shell type(s) {bad} not "
        "supported by this function (only S, P, SP, D, F, and G -- pure or "
        "Cartesian -- shells are implemented). Rendering would be silently wrong "
        "for any basis function belonging to these shells, so no surface is "
        "produced. A common cause is an h-polarized or larger basis set (rare) "
        "or a non-standard shell-type code."
    )


# =============================================================================
# .fchk raw-section parsing (self-contained -- does NOT reuse G16parser's
# fchk_read, mirroring the MATLAB original's own deliberate independence
# from G09/G16_fchk_read)
# =============================================================================

_HEADER_RE = re.compile(r'^(.{1,43}?)\s{1,3}(I|R|C)\s+(N=\s*\d+|[-\d.E+]+)')


def _read_fchk_lines(filename):
    with open(filename, "r", encoding="utf-8", errors="replace") as f:
        return f.read().split("\n")


def _index_fchk_sections(lines):
    """Returns a list of (name, dtype, nvals, line_index) for every
    scalar/array section header found in lines.
    """
    sections = []
    for k, ln in enumerate(lines):
        if len(ln) < 45:
            continue
        m = _HEADER_RE.match(ln)
        if not m:
            continue
        name = m.group(1).strip()
        dtype = m.group(2)
        rest = m.group(3).strip()
        if rest.startswith("N="):
            nvals = int(rest[2:].strip())
        else:
            nvals = 0
        sections.append((name, dtype, nvals, k))
    return sections


def _read_fchk_section(lines, sections, keyword):
    """Finds keyword (case-insensitive, last match wins, matching the
    MATLAB original) among sections and reads its value(s): a bare float
    for a scalar section, or a numpy array for an array section. Returns
    None if not found.
    """
    matches = [s for s in sections if s[0].lower() == keyword.lower()]
    if not matches:
        return None
    name, dtype, nvals, k0 = matches[-1]
    N = len(lines)

    if nvals == 0:
        m = re.search(r'(I|R|C)\s+([-\d.E+]+)\s*$', lines[k0])
        if not m:
            return None
        return float(m.group(2))

    vals = []
    k2 = k0 + 1
    while len(vals) < nvals and k2 < N:
        ln2 = lines[k2].strip()
        if not ln2:
            k2 += 1
            continue
        if len(ln2) > 44 and _HEADER_RE.match(ln2):
            break
        vals.extend(float(t) for t in ln2.split())
        k2 += 1
    return np.array(vals[:nvals])


def read_aobasis_from_fchk(filename):
    """Minimal, standalone .fchk section parser that reads only the raw
    atomic-orbital basis-set sections (shell types, primitive exponents,
    contraction coefficients, shell centres, ...) needed to evaluate a
    molecular orbital or the electron density on a real-space grid.

    Returns a dict with keys shell_types, n_prim_per_shell, shell_to_atom
    (all int arrays), prim_exponents, contraction_coeff,
    sp_contraction_coeff (None if no SP shells), shell_coords_bohr
    ([Nshell, 3], None if the section was absent/malformed).
    """
    lines = _read_fchk_lines(filename)
    sections = _index_fchk_sections(lines)

    shell_types = _read_fchk_section(lines, sections, "Shell types")
    n_prim_per_shell = _read_fchk_section(lines, sections, "Number of primitives per shell")
    shell_to_atom = _read_fchk_section(lines, sections, "Shell to atom map")
    prim_exponents = _read_fchk_section(lines, sections, "Primitive exponents")
    contraction_coeff = _read_fchk_section(lines, sections, "Contraction coefficients")
    sp_contraction_coeff = _read_fchk_section(lines, sections, "P(S=P) Contraction coefficients")
    shell_coords_raw = _read_fchk_section(lines, sections, "Coordinates of each shell")

    aobasis = {
        "shell_types": np.round(shell_types).astype(int) if shell_types is not None else np.array([], dtype=int),
        "n_prim_per_shell": np.round(n_prim_per_shell).astype(int) if n_prim_per_shell is not None else np.array([], dtype=int),
        "shell_to_atom": np.round(shell_to_atom).astype(int) if shell_to_atom is not None else np.array([], dtype=int),
        "prim_exponents": prim_exponents if prim_exponents is not None else np.array([]),
        "contraction_coeff": contraction_coeff if contraction_coeff is not None else np.array([]),
        "sp_contraction_coeff": sp_contraction_coeff,
    }
    nshell = aobasis["shell_types"].size
    if shell_coords_raw is not None and shell_coords_raw.size == 3 * nshell:
        aobasis["shell_coords_bohr"] = shell_coords_raw.reshape(nshell, 3)
    else:
        aobasis["shell_coords_bohr"] = None
    return aobasis


def read_fchk_scalar_or_section(filename, keyword):
    """Public convenience: read one named .fchk section (scalar float or
    array) directly, e.g. "Nuclear charges" or "Beta MO coefficients" --
    used by the ESP/spin-density code paths for sections outside the
    basis-set scope of read_aobasis_from_fchk.
    """
    lines = _read_fchk_lines(filename)
    sections = _index_fchk_sections(lines)
    return _read_fchk_section(lines, sections, keyword)


# =============================================================================
# Gaussian .cube file I/O
# =============================================================================

def write_cube_file(filename, title_line, comment_line, atomic_numbers, xyz_bohr,
                     gx, gy, gz, V, orbital_indices=None):
    """Writes a scalar volumetric grid to a Gaussian-format .cube file.

    atomic_numbers : [Natoms] true atomic numbers (not ECP-reduced charges)
    xyz_bohr       : [Natoms, 3]
    gx, gy, gz     : 1D axis vectors (Bohr), uniformly spaced
    V              : meshgrid convention (numpy 'ij'... see note below):
                     V has shape (len(gy), len(gx), len(gz)) -- i.e.
                     V[iy, ix, iz] -- the SAME convention produced by
                     np.meshgrid(gx, gy, gz) (default 'xy' indexing) and
                     used internally by every g16_draw_*_surface caller.
    orbital_indices - optional 1-based MO index list; if given, Natoms is
                     written as NEGATIVE and an extra orbital-index line
                     follows the atom list (Gaussian MO-cube convention).

    Value ordering: a cube file's flat data is X slowest, Z fastest.
    Given V(iy,ix,iz), the flat order is V.transpose(2,0,1).reshape(-1) --
    the same recipe as the validated MATLAB original
    (permute(V,[3 1 2])), just 0-based.
    """
    atomic_numbers = np.asarray(atomic_numbers)
    Natoms = atomic_numbers.size
    Nx, Ny, Nz = len(gx), len(gy), len(gz)
    if V.shape != (Ny, Nx, Nz):
        raise ValueError(
            f"write_cube_file: V has shape {V.shape}, expected ({Ny}, {Nx}, {Nz}) "
            "= (len(gy), len(gx), len(gz)) (meshgrid convention)."
        )

    dx = gx[1] - gx[0] if Nx > 1 else 1.0
    dy = gy[1] - gy[0] if Ny > 1 else 1.0
    dz = gz[1] - gz[0] if Nz > 1 else 1.0
    origin_bohr = (gx[0], gy[0], gz[0])

    is_mo_cube = orbital_indices is not None and len(orbital_indices) > 0
    natoms_field = -Natoms if is_mo_cube else Natoms

    with open(filename, "w") as f:
        f.write(f"{title_line}\n")
        f.write(f"{comment_line}\n")
        f.write(f"{natoms_field:5d}{origin_bohr[0]:12.6f}{origin_bohr[1]:12.6f}{origin_bohr[2]:12.6f}\n")
        f.write(f"{Nx:5d}{dx:12.6f}{0.0:12.6f}{0.0:12.6f}\n")
        f.write(f"{Ny:5d}{0.0:12.6f}{dy:12.6f}{0.0:12.6f}\n")
        f.write(f"{Nz:5d}{0.0:12.6f}{0.0:12.6f}{dz:12.6f}\n")

        for A in range(Natoms):
            z = atomic_numbers[A]
            f.write(f"{int(z):5d}{float(z):12.6f}{xyz_bohr[A,0]:12.6f}{xyz_bohr[A,1]:12.6f}{xyz_bohr[A,2]:12.6f}\n")

        if is_mo_cube:
            f.write(f"{len(orbital_indices):5d}")
            for oi in orbital_indices:
                f.write(f"{int(oi):5d}")
            f.write("\n")

        # Flatten to Gaussian cube order (X slowest, Y next, Z fastest).
        # NOTE: MATLAB is column-major, numpy is row-major -- the MATLAB
        # original's permute(V,[3 1 2]) + reshape(...,Nz,[]) recipe is
        # NOT a literal transliteration target here; re-derived natively
        # instead. V has shape (Ny,Nx,Nz); P=transpose(V,(2,1,0)) has
        # shape (Nz,Nx,Ny) with P[iz,ix,iy]=V[iy,ix,iz]. A row-major
        # reshape(Nz, Nx*Ny) then combines the trailing (Nx,Ny) axes with
        # Ny (iy) fastest and Nx (ix) slower -- column index c=ix*Ny+iy,
        # i.e. ix outer/slowest, iy inner -- exactly the required z-run
        # enumeration order. Verified against an explicit triple nested
        # loop (ix outer, iy, iz inner) on a small test grid, exact match.
        P = np.transpose(V, (2, 1, 0))
        zruns = P.reshape(Nz, Nx * Ny)

        # A real cube file force-breaks the line at the END of every
        # z-run, not just every 6 values continuously across the whole
        # flat stream (see g_write_cube_file.m's header for the
        # GaussView-crash story this fixes).
        for r in range(zruns.shape[1]):
            col = zruns[:, r]
            for k0 in range(0, Nz, 6):
                chunk = col[k0:k0 + 6]
                f.write("".join(f"{v:13.5E}" for v in chunk))
                f.write("\n")


def read_cube_file(filename):
    """Reads a Gaussian-format .cube file (the inverse of write_cube_file,
    compatible with real cubegen/GaussView output too).

    Returns a dict with keys: atomic_numbers [Natoms], xyz_bohr
    [Natoms,3], gx, gy, gz (1D Bohr axis vectors), V (meshgrid
    convention, V[iy,ix,iz]), orbital_indices (None for a scalar-field
    cube), title_line, comment_line.
    """
    with open(filename, "r", encoding="utf-8", errors="replace") as f:
        lines = f.read().split("\n")

    title_line = lines[0]
    comment_line = lines[1]

    l3 = [float(t) for t in lines[2].split()]
    natoms_field = round(l3[0])
    origin_bohr = np.array(l3[1:4])
    is_mo_cube = natoms_field < 0
    natoms = abs(natoms_field)

    N = [0, 0, 0]
    step = np.zeros((3, 3))
    for k in range(3):
        lk = [float(t) for t in lines[3 + k].split()]
        N[k] = round(lk[0])
        step[k, :] = lk[1:4]
    off_diag = step - np.diag(np.diag(step))
    if np.any(np.abs(off_diag) > 1e-9):
        raise ValueError(
            f"read_cube_file: {filename} has a non-axis-aligned grid (off-diagonal "
            "step vectors) -- not supported by this toolbox's grid-based tools."
        )

    row = 6
    atomic_numbers = np.zeros(natoms, dtype=int)
    xyz_bohr = np.zeros((natoms, 3))
    for a in range(natoms):
        la = [float(t) for t in lines[row].split()]
        atomic_numbers[a] = round(la[0])
        xyz_bohr[a, :] = la[2:5]
        row += 1

    orbital_indices = None
    if is_mo_cube:
        ln = [float(t) for t in lines[row].split()]
        row += 1
        n_orb = round(ln[0])
        orbital_indices = np.round(ln[1:1 + n_orb]).astype(int)

    # Remaining lines are the flat data stream (whitespace/newline
    # agnostic, unlike MATLAB's token-based fscanf equivalent here).
    rest_text = "\n".join(lines[row:])
    vals = np.array([float(t) for t in rest_text.split()])

    Nx, Ny, Nz = N
    if vals.size != Nx * Ny * Nz:
        raise ValueError(
            f"read_cube_file: {filename}: read {vals.size} data values, expected "
            f"{Nx*Ny*Nz} ({Nx} x {Ny} x {Nz})."
        )

    # Flat order is X slowest, Y next, Z fastest -- i.e. exactly the
    # nested-loop order (ix outer, iy, iz inner) a row-major/C-order
    # reshape(Nx,Ny,Nz) reconstructs directly (last axis fastest, first
    # axis slowest -- no permute needed, unlike the MATLAB original's
    # column-major-specific two-step reshape+permute).
    Vxyz = vals.reshape(Nx, Ny, Nz)   # Vxyz[ix,iy,iz]
    V = np.transpose(Vxyz, (1, 0, 2))  # V[iy,ix,iz], meshgrid convention

    gx = origin_bohr[0] + np.arange(Nx) * step[0, 0]
    gy = origin_bohr[1] + np.arange(Ny) * step[1, 1]
    gz = origin_bohr[2] + np.arange(Nz) * step[2, 2]

    return {
        "atomic_numbers": atomic_numbers,
        "xyz_bohr": xyz_bohr,
        "gx": gx, "gy": gy, "gz": gz,
        "V": V,
        "orbital_indices": orbital_indices,
        "title_line": title_line,
        "comment_line": comment_line,
    }


def cube_grid_params(cubefile):
    """Recovers the grid spacing and padding used to build a .cube file,
    so a companion cube (e.g. for g16_draw_cube_surface's color_by,
    which requires two cubes to share an identical grid) can be
    regenerated with a matching grid_spacing/padding (or
    cube_spacing/cube_padding for g16_draw_esp_surface's own save_cube
    grid) without having to note them down manually.

    Works on any axis-aligned, uniformly-spaced .cube file containing
    atoms -- one written by this package's write_cube_file (via any
    g16_draw_*_surface save_cube option), or a real cubegen/GaussView
    cube.

    Returns a dict with keys:
        grid_spacing      - float, Bohr (the cube's step size; raises
                             ValueError if dx/dy/dz differ by more than
                             1e-6 Bohr, i.e. the cube is not isotropic)
        padding            - float, Bohr: the mean of the six per-face
                             paddings below -- the value to pass back as
                             padding/cube_padding for a matching cube
        padding_per_face   - (2,3) ndarray, Bohr: [low;high] x [x,y,z],
                             the raw per-face (grid boundary minus atom
                             bounding box) distances padding was averaged
                             from. A warning is printed if these differ
                             by more than one grid_spacing, since this
                             package's own grids always use the SAME
                             padding on every side -- a wide spread here
                             means the cube was not built that way, and
                             padding alone will not exactly reproduce it.
        npoints            - (Nx, Ny, Nz)
        origin_bohr        - (3,) ndarray

    CAUTION -- direction matters: ESP is far more expensive to evaluate
    per grid point than density or an MO (McMurchie-Davidson integrals,
    or a Coulomb-grid-sum fallback for basis sets with shells the
    analytic method does not support), and g16_draw_esp_surface's own
    save_cube path evaluates ESP at every point of the saved grid, not
    just at isosurface vertices, capped at 2e6 points so it errors out
    instead of silently running for hours. A grid recovered from a
    DENSITY cube (typically 0.10-0.15 Bohr spacing, fine enough for a
    nice-looking isosurface) is very often too fine to reuse for an ESP
    save_cube this way. Prefer the opposite order: pick a shared
    grid_spacing/padding the ESP evaluation can afford (e.g. around
    g16_draw_esp_surface's own defaults, cube_spacing=0.30/
    cube_padding=4.0) for the ESP cube first, then match the (much
    cheaper to re-evaluate) density cube to IT:

        g16_draw_esp_surface(data, save_cube='esp.cube',
            cube_spacing=0.30, cube_padding=4.0)
        p = cube_grid_params('esp.cube')
        g16_draw_density_surface(data, save_cube='density.cube',
            grid_spacing=p['grid_spacing'], padding=p['padding'])
        g16_draw_cube_surface('density.cube', color_by='esp.cube')
    """
    d = read_cube_file(cubefile)
    gx, gy, gz = d["gx"], d["gy"], d["gz"]
    xyz_bohr = d["xyz_bohr"]

    dx = gx[1] - gx[0]
    dy = gy[1] - gy[0]
    dz = gz[1] - gz[0]
    if max(abs(dx - dx), abs(dy - dx), abs(dz - dx)) > 1e-6:
        raise ValueError(
            f"cube_grid_params: {cubefile} has different spacing along x/y/z "
            f"({dx:.6g}/{dy:.6g}/{dz:.6g} Bohr) -- not an isotropic grid as "
            "written by this package."
        )
    spacing = dx

    if xyz_bohr.shape[0] == 0:
        raise ValueError(
            f"cube_grid_params: {cubefile} has no atoms -- cannot recover "
            "padding (only grid_spacing)."
        )

    atom_lo = xyz_bohr.min(axis=0)
    atom_hi = xyz_bohr.max(axis=0)
    grid_lo = np.array([gx[0], gy[0], gz[0]])
    grid_hi = np.array([gx[-1], gy[-1], gz[-1]])

    # pad_lo is EXACT by construction: every g16_draw_*_surface/
    # write_cube_file caller builds its grid as grid_lo = atom_lo -
    # padding, so this recovers the true original value. pad_hi is only
    # an underestimate (by up to one grid step), since the numpy.arange
    # grid stops at the last point at or below grid_hi -- expected, not
    # an error, and NOT used for the returned padding (using it, or its
    # mean with pad_lo, would fail to reproduce the same grid_lo when
    # fed back in, breaking g16_draw_cube_surface's color_by).
    pad_lo = atom_lo - grid_lo
    pad_hi = grid_hi - atom_hi
    pad_all = np.vstack([pad_lo, pad_hi])   # (2,3), kept for diagnostics only

    padding = float(pad_lo.mean())
    if pad_lo.max() - pad_lo.min() > spacing:
        print(
            f"cube_grid_params: WARNING {cubefile}: the low-side padding varies "
            f"by more than one grid step across x/y/z ({pad_lo.min():.3g} to "
            f"{pad_lo.max():.3g} Bohr) -- this cube was likely not built with a "
            f"single uniform padding value; the returned padding ({padding:.3g} "
            "Bohr) may not exactly reproduce this grid."
        )

    return {
        "grid_spacing": spacing,
        "padding": padding,
        "padding_per_face": pad_all,
        "npoints": (len(gx), len(gy), len(gz)),
        "origin_bohr": grid_lo,
    }
