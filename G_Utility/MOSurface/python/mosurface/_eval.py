"""Real-space Gaussian-type-orbital evaluation engine: molecular orbitals,
electron density, and the analytic (McMurchie-Davidson) electrostatic
potential -- the numerical core shared by all four g16_draw_*_surface
functions. Direct port of g_eval_mo_on_grid.m / g_eval_density_on_grid.m /
g_eval_esp_analytic.m; see those files' headers for the full derivation
and validation notes (normalization conventions, shell orderings,
McMurchie-Davidson references).
"""

import numpy as np

from ._common import boys_function


# =============================================================================
# Odd double factorial and normalization helpers
# =============================================================================

def _ddf(p):
    """Odd double factorial (2p-1)!!, ddf(0) = 1 by convention. Accepts a
    non-integer p (only ever arising, harmlessly, as an intermediate half-
    sum in a self-overlap computation whose paired-monomial-degree parity
    always works out even for every actual use in this module): MATLAB's
    original used prod(1:2:(2p-1)), whose colon range naturally floors to
    floor(p) for a non-integer p, replicated here explicitly.
    """
    k = int(np.floor(p))
    if k <= 0:
        return 1
    w = 1
    for x in range(1, 2 * k, 2):
        w *= x
    return w


def _norm_const_ijk(alpha, i, j, k):
    """Individual primitive normalization for the Cartesian monomial
    x^i y^j z^k: N = (2a/pi)^(3/4) sqrt((4a)^L / (ddf(i)ddf(j)ddf(k))).
    Must be evaluated per monomial, not once per shell from L alone (see
    g_eval_esp_analytic.m's header for the bug this distinction fixes).
    """
    L = i + j + k
    return (2 * alpha / np.pi) ** 0.75 * np.sqrt((4 * alpha) ** L / (_ddf(i) * _ddf(j) * _ddf(k)))


# =============================================================================
# Pure D/F/G real solid harmonics (Ribaldone & Desmarais 2024, Table I)
# =============================================================================

def _f_pure_terms():
    return [
        [(1, 0, 0, 3), (-1.5, 2, 0, 1), (-1.5, 0, 2, 1)],       # m=0
        [(6, 1, 0, 2), (-1.5, 1, 2, 0), (-1.5, 3, 0, 0)],       # m=+1
        [(6, 0, 1, 2), (-1.5, 0, 3, 0), (-1.5, 2, 1, 0)],       # m=-1
        [(-15, 0, 2, 1), (15, 2, 0, 1)],                        # m=+2
        [(30, 1, 1, 1)],                                        # m=-2
        [(-45, 1, 2, 0), (15, 3, 0, 0)],                        # m=+3
        [(-15, 0, 3, 0), (45, 2, 1, 0)],                        # m=-3
    ]


def _g_pure_terms():
    return [
        [(1, 0, 0, 4), (-3, 0, 2, 2), (0.375, 0, 4, 0), (-3, 2, 0, 2), (0.75, 2, 2, 0), (0.375, 4, 0, 0)],  # m=0
        [(10, 1, 0, 3), (-7.5, 1, 2, 1), (-7.5, 3, 0, 1)],       # m=+1
        [(10, 0, 1, 3), (-7.5, 0, 3, 1), (-7.5, 2, 1, 1)],       # m=-1
        [(-45, 0, 2, 2), (7.5, 0, 4, 0), (45, 2, 0, 2), (-7.5, 4, 0, 0)],  # m=+2
        [(90, 1, 1, 2), (-15, 1, 3, 0), (-15, 3, 1, 0)],         # m=-2
        [(-315, 1, 2, 1), (105, 3, 0, 1)],                       # m=+3
        [(-105, 0, 3, 1), (315, 2, 1, 1)],                       # m=-3
        [(105, 0, 4, 0), (-630, 2, 2, 0), (105, 4, 0, 0)],       # m=+4
        [(-420, 1, 3, 0), (420, 3, 1, 0)],                       # m=-4
    ]


def _pure_d_components():
    """5 pure-D basis functions (order z2,xz,yz,x2-y2,xy), each a linear
    combination of individually-normalized Cartesian D monomials, with
    the overall normalization K folded into the returned coefficients
    (used by the analytic ESP evaluator; see g_eval_esp_analytic.m's
    header for the individually-normalized-overlap derivation this
    needs -- NOT a single shared (2L-1)!! divisor).
    """
    terms = [
        [(1, 0, 0, 2), (-0.5, 2, 0, 0), (-0.5, 0, 2, 0)],           # m=0 (z2)
        [(1, 1, 0, 1)],                                              # m=+1 (xz)
        [(1, 0, 1, 1)],                                              # m=-1 (yz)
        [(0.5 * np.sqrt(3), 2, 0, 0), (-0.5 * np.sqrt(3), 0, 2, 0)],  # m=+2 (x2-y2)
        [(1, 1, 1, 0)],                                               # m=-2 (xy)
    ]
    components = []
    for term in terms:
        P = np.array([t[0] for t in term])
        mono = np.array([t[1:] for t in term])
        nterm = len(term)
        D = np.array([_ddf(mono[i, 0]) * _ddf(mono[i, 1]) * _ddf(mono[i, 2]) for i in range(nterm)])
        S = 0.0
        for i in range(nterm):
            for j in range(nterm):
                raw = (_ddf((mono[i, 0] + mono[j, 0]) / 2) *
                       _ddf((mono[i, 1] + mono[j, 1]) / 2) *
                       _ddf((mono[i, 2] + mono[j, 2]) / 2))
                S += P[i] * P[j] * raw / np.sqrt(D[i] * D[j])
        K = 1.0 / np.sqrt(S)
        components.append({"mono": mono, "coef": K * P})
    return components


def _pure_harmonic_value(terms, dx, dy, dz, radial_shared, L):
    """Evaluates one normalized real solid harmonic pure D/F/G component
    from its unnormalized polynomial terms=[(P,t,u,v), ...], combined
    with radial_shared (normalized for the "diagonal" x^L component). The
    overall normalization K is computed at runtime from the exact
    combinatorial self-overlap of same-center Cartesian Gaussians.
    """
    P = np.array([t[0] for t in terms])
    tt = np.array([t[1] for t in terms])
    uu = np.array([t[2] for t in terms])
    vv = np.array([t[3] for t in terms])
    n = len(terms)
    S = 0.0
    for k in range(n):
        for l in range(n):
            S += P[k] * P[l] * _ddf((tt[k] + tt[l]) // 2) * _ddf((uu[k] + uu[l]) // 2) * _ddf((vv[k] + vv[l]) // 2)
    K = np.sqrt(_ddf(L) / S)
    poly = np.zeros_like(dx)
    for k in range(n):
        poly = poly + P[k] * (dx ** tt[k]) * (dy ** uu[k]) * (dz ** vv[k])
    return K * poly * radial_shared


def _cart_component_value(t, u, v, dx, dy, dz, radial_shared, L):
    """Value of one Cartesian component (t,u,v), t+u+v=L, of a Cartesian
    D/F/G shell, correcting radial_shared (normalized for the "diagonal"
    component) to this component's own individual normalization.
    """
    corr = np.sqrt(_ddf(L) / (_ddf(t) * _ddf(u) * _ddf(v)))
    return corr * (dx ** t) * (dy ** u) * (dz ** v) * radial_shared


def _local_radial(alphas, coeffs, r2, L):
    """Contracted radial part shared by all Cartesian components of
    angular momentum L within one shell (using the "diagonal" primitive
    normalization).
    """
    dfact = _ddf(L)
    Nprim = (2 * alphas / np.pi) ** 0.75 * np.sqrt((4 * alphas) ** L / dfact)
    radial = np.zeros_like(r2)
    for pk in range(len(alphas)):
        radial = radial + coeffs[pk] * Nprim[pk] * np.exp(-alphas[pk] * r2)
    return radial


# =============================================================================
# eval_mo_on_grid / eval_density_on_grid
# =============================================================================

def eval_mo_on_grid(aobasis, mo_coeff, grid_pts_bohr):
    """Evaluates one or more MOs' real-space wavefunctions on a grid of
    points [Ngrid, 3] (Bohr), accumulating shell-by-shell directly into
    an [Ngrid, K] array -- memory stays O(Ngrid*K) regardless of basis
    set size.

    mo_coeff : [Nbasis, K] MO coefficient columns (K=1 for a single MO).
    Returns psi: [Ngrid, K].
    """
    grid_pts_bohr = np.asarray(grid_pts_bohr)
    mo_coeff = np.atleast_2d(mo_coeff)
    if mo_coeff.shape[0] == 1 and mo_coeff.shape[1] != 1:
        mo_coeff = mo_coeff.T  # allow a plain 1D vector to mean a single MO column
    Ngrid = grid_pts_bohr.shape[0]
    K = mo_coeff.shape[1]
    shell_types = aobasis["shell_types"]
    n_prim_per_shell = aobasis["n_prim_per_shell"]
    Nshell = len(shell_types)
    psi = np.zeros((Ngrid, K))

    prim_offset = 0
    col = 0
    for s in range(Nshell):
        t = shell_types[s]
        nprim = n_prim_per_shell[s]
        idx_p = slice(prim_offset, prim_offset + nprim)
        alphas = aobasis["prim_exponents"][idx_p]
        cS = aobasis["contraction_coeff"][idx_p]
        center = aobasis["shell_coords_bohr"][s, :]

        dx = grid_pts_bohr[:, 0] - center[0]
        dy = grid_pts_bohr[:, 1] - center[1]
        dz = grid_pts_bohr[:, 2] - center[2]
        r2 = dx**2 + dy**2 + dz**2

        if t == 0:  # S
            radial = _local_radial(alphas, cS, r2, 0)
            psi += np.outer(radial, mo_coeff[col, :])
            col += 1

        elif t == 1:  # P: X,Y,Z
            radial = _local_radial(alphas, cS, r2, 1)
            psi += (np.outer(radial * dx, mo_coeff[col, :])
                    + np.outer(radial * dy, mo_coeff[col + 1, :])
                    + np.outer(radial * dz, mo_coeff[col + 2, :]))
            col += 3

        elif t == -1:  # SP: S, X, Y, Z
            cP = aobasis["sp_contraction_coeff"][idx_p]
            radialS = _local_radial(alphas, cS, r2, 0)
            radialP = _local_radial(alphas, cP, r2, 1)
            psi += np.outer(radialS, mo_coeff[col, :])
            psi += (np.outer(radialP * dx, mo_coeff[col + 1, :])
                    + np.outer(radialP * dy, mo_coeff[col + 2, :])
                    + np.outer(radialP * dz, mo_coeff[col + 3, :]))
            col += 4

        elif t == -2:  # pure D: z2, xz, yz, x2-y2, xy
            radialD = _local_radial(alphas, cS, r2, 2)
            Dxx = dx * dx * radialD
            Dyy = dy * dy * radialD
            Dzz = dz * dz * radialD
            Dxy = dx * dy * radialD * np.sqrt(3)
            Dxz = dx * dz * radialD * np.sqrt(3)
            Dyz = dy * dz * radialD * np.sqrt(3)
            d_z2 = Dzz - 0.5 * Dxx - 0.5 * Dyy
            d_x2y2 = (np.sqrt(3) / 2) * (Dxx - Dyy)
            psi += (np.outer(d_z2, mo_coeff[col, :]) + np.outer(Dxz, mo_coeff[col + 1, :])
                    + np.outer(Dyz, mo_coeff[col + 2, :]) + np.outer(d_x2y2, mo_coeff[col + 3, :])
                    + np.outer(Dxy, mo_coeff[col + 4, :]))
            col += 5

        elif t == 2:  # Cartesian D: xx,yy,zz,xy,xz,yz
            radialD = _local_radial(alphas, cS, r2, 2)
            psi += (np.outer(dx * dx * radialD, mo_coeff[col, :])
                    + np.outer(dy * dy * radialD, mo_coeff[col + 1, :])
                    + np.outer(dz * dz * radialD, mo_coeff[col + 2, :])
                    + np.outer(np.sqrt(3) * dx * dy * radialD, mo_coeff[col + 3, :])
                    + np.outer(np.sqrt(3) * dx * dz * radialD, mo_coeff[col + 4, :])
                    + np.outer(np.sqrt(3) * dy * dz * radialD, mo_coeff[col + 5, :]))
            col += 6

        elif t == -3:  # pure F (7 components)
            radialF = _local_radial(alphas, cS, r2, 3)
            for c, terms in enumerate(_f_pure_terms()):
                val = _pure_harmonic_value(terms, dx, dy, dz, radialF, 3)
                psi += np.outer(val, mo_coeff[col + c, :])
            col += 7

        elif t == 3:  # Cartesian F (10), Gaussian native order
            radialF = _local_radial(alphas, cS, r2, 3)
            cart_f = [(3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0),
                      (2, 0, 1), (1, 0, 2), (0, 1, 2), (0, 2, 1), (1, 1, 1)]
            for c, (i_, j_, k_) in enumerate(cart_f):
                val = _cart_component_value(i_, j_, k_, dx, dy, dz, radialF, 3)
                psi += np.outer(val, mo_coeff[col + c, :])
            col += 10

        elif t == -4:  # pure G (9 components)
            radialG = _local_radial(alphas, cS, r2, 4)
            for c, terms in enumerate(_g_pure_terms()):
                val = _pure_harmonic_value(terms, dx, dy, dz, radialG, 4)
                psi += np.outer(val, mo_coeff[col + c, :])
            col += 9

        elif t == 4:  # Cartesian G (15), Gaussian native (reversed-lex) order
            radialG = _local_radial(alphas, cS, r2, 4)
            cart_g = [(0, 0, 4), (0, 1, 3), (0, 2, 2), (0, 3, 1), (0, 4, 0),
                      (1, 0, 3), (1, 1, 2), (1, 2, 1), (1, 3, 0),
                      (2, 0, 2), (2, 1, 1), (2, 2, 0), (3, 0, 1), (3, 1, 0), (4, 0, 0)]
            for c, (i_, j_, k_) in enumerate(cart_g):
                val = _cart_component_value(i_, j_, k_, dx, dy, dz, radialG, 4)
                psi += np.outer(val, mo_coeff[col + c, :])
            col += 15

        prim_offset += nprim

    return psi


def eval_density_on_grid(aobasis, occ_coeff, grid_pts_bohr, occ_factor=2):
    """rho(r) = occ_factor * sum_i psi_i(r)^2 over the occupied orbitals
    in occ_coeff (an [Nbasis, Nocc] matrix), evaluated in a single pass
    through the basis shells. Returns rho: [Ngrid].
    """
    psi_occ = eval_mo_on_grid(aobasis, occ_coeff, grid_pts_bohr)
    return occ_factor * np.sum(psi_occ**2, axis=1)


# =============================================================================
# Analytic (McMurchie-Davidson) electrostatic potential
# =============================================================================

def _s_components():
    return [{"mono": np.array([[0, 0, 0]]), "coef": np.array([1.0])}]


def _p_components():
    return [
        {"mono": np.array([[1, 0, 0]]), "coef": np.array([1.0])},
        {"mono": np.array([[0, 1, 0]]), "coef": np.array([1.0])},
        {"mono": np.array([[0, 0, 1]]), "coef": np.array([1.0])},
    ]


def _hermite_E_1d(alpha, beta, AX, BX, La, Lb):
    """McMurchie-Davidson Hermite expansion coefficients E_t^{ij}, one
    Cartesian direction, i=0..La, j=0..Lb, t=0..i+j. 0-based array E of
    shape (La+1, Lb+1, La+Lb+1) -- E[i,j,t] (no MATLAB +1 offsets needed,
    numpy is already 0-based).
    """
    p = alpha + beta
    mu = alpha * beta / p
    P = (alpha * AX + beta * BX) / p
    XPA = P - AX
    XPB = P - BX
    XAB = AX - BX
    Tmax = La + Lb
    E = np.zeros((La + 1, Lb + 1, Tmax + 1))
    E[0, 0, 0] = np.exp(-mu * XAB**2)

    for j in range(1, Lb + 1):
        ij_prev = j - 1
        for t in range(0, j + 1):
            term1 = E[0, j - 1, t - 1] if (t >= 1 and (t - 1) <= ij_prev) else 0.0
            term2 = E[0, j - 1, t] if (t <= ij_prev) else 0.0
            term3 = E[0, j - 1, t + 1] if ((t + 1) <= ij_prev) else 0.0
            E[0, j, t] = term1 / (2 * p) + XPB * term2 + (t + 1) * term3

    for i in range(1, La + 1):
        for j in range(0, Lb + 1):
            ij_prev = (i - 1) + j
            for t in range(0, i + j + 1):
                term1 = E[i - 1, j, t - 1] if (t >= 1 and (t - 1) <= ij_prev) else 0.0
                term2 = E[i - 1, j, t] if (t <= ij_prev) else 0.0
                term3 = E[i - 1, j, t + 1] if ((t + 1) <= ij_prev) else 0.0
                E[i, j, t] = term1 / (2 * p) + XPA * term2 + (t + 1) * term3

    return E


def _hermite_R_full(p, PC, Ltot):
    """McMurchie-Davidson Hermite Coulomb integrals R_tuv (n=0), for all
    t+u+v<=Ltot, vectorized over field points. PC: [Nfield, 3] (P - C for
    each field point C). Returns R: [Ltot+1, Ltot+1, Ltot+1, Nfield].
    """
    Nfield = PC.shape[0]
    Rfull = np.zeros((Ltot + 1, Ltot + 1, Ltot + 1, Ltot + 1, Nfield))

    x = p * np.sum(PC**2, axis=1)
    Fn = boys_function(x, Ltot)   # [Nfield, Ltot+1]
    for n in range(Ltot + 1):
        Rfull[0, 0, 0, n, :] = ((-2 * p) ** n) * Fn[:, n]

    for D in range(1, Ltot + 1):
        for t in range(0, D + 1):
            for u in range(0, D - t + 1):
                v = D - t - u
                if t > 0:
                    for n in range(0, Ltot - D + 1):
                        term = PC[:, 0] * Rfull[t - 1, u, v, n + 1, :]
                        if t >= 2:
                            term = term + (t - 1) * Rfull[t - 2, u, v, n + 1, :]
                        Rfull[t, u, v, n, :] = term
                elif u > 0:
                    for n in range(0, Ltot - D + 1):
                        term = PC[:, 1] * Rfull[0, u - 1, v, n + 1, :]
                        if u >= 2:
                            term = term + (u - 1) * Rfull[0, u - 2, v, n + 1, :]
                        Rfull[0, u, v, n, :] = term
                else:
                    for n in range(0, Ltot - D + 1):
                        term = PC[:, 2] * Rfull[0, 0, v - 1, n + 1, :]
                        if v >= 2:
                            term = term + (v - 1) * Rfull[0, 0, v - 2, n + 1, :]
                        Rfull[0, 0, v, n, :] = term

    return Rfull[:, :, :, 0, :]


class _Subshell:
    __slots__ = ("L", "alphas", "coeffs", "center", "components", "col0")

    def __init__(self, L, alphas, coeffs, center, components, col0):
        self.L = L
        self.alphas = alphas
        self.coeffs = coeffs
        self.center = center
        self.components = components
        self.col0 = col0


def eval_esp_analytic(aobasis, P_density, field_pts_bohr, verbose=False):
    """Exact analytic electronic-Coulomb-potential evaluation
        V_elec(r) = sum_munu P_munu <phi_mu | 1/|r-r'| | phi_nu>
    at each field point r, via the McMurchie-Davidson Hermite-Gaussian
    scheme. Only S, P, SP, and pure D shells are supported (matching
    g_eval_esp_analytic.m's scope); raises ValueError otherwise.

    P_density : [Nbasis, Nbasis] density matrix.
    field_pts_bohr : [Nfield, 3].
    Returns V_elec: [Nfield].
    """
    field_pts_bohr = np.asarray(field_pts_bohr)
    Nfield = field_pts_bohr.shape[0]
    V_elec = np.zeros(Nfield)
    shell_types = aobasis["shell_types"]
    n_prim_per_shell = aobasis["n_prim_per_shell"]
    Nshell = len(shell_types)

    subshells = []
    prim_offset = 0
    col = 0
    for s in range(Nshell):
        t = shell_types[s]
        nprim = n_prim_per_shell[s]
        idx_p = slice(prim_offset, prim_offset + nprim)
        alphas = aobasis["prim_exponents"][idx_p]
        cS = aobasis["contraction_coeff"][idx_p]
        center = aobasis["shell_coords_bohr"][s, :]

        if t == 0:
            subshells.append(_Subshell(0, alphas, cS, center, _s_components(), col))
            col += 1
        elif t == 1:
            subshells.append(_Subshell(1, alphas, cS, center, _p_components(), col))
            col += 3
        elif t == -1:
            cP = aobasis["sp_contraction_coeff"][idx_p]
            subshells.append(_Subshell(0, alphas, cS, center, _s_components(), col))
            col += 1
            subshells.append(_Subshell(1, alphas, cP, center, _p_components(), col))
            col += 3
        elif t == -2:
            subshells.append(_Subshell(2, alphas, cS, center, _pure_d_components(), col))
            col += 5
        else:
            raise ValueError(
                f"eval_esp_analytic: only supports S, P, SP, and pure D shells "
                f"(shell type {t} at shell {s} is not one of these). Cartesian D and "
                "F/G shells are not yet implemented for the analytic ESP evaluator."
            )
        prim_offset += nprim

    if verbose:
        print(f"eval_esp_analytic: {len(subshells)} subshells, {col} basis functions, {Nfield} field points")

    Nsub = len(subshells)
    for iA in range(Nsub):
        A = subshells[iA]
        ncompA = len(A.components)
        for iB in range(Nsub):
            B = subshells[iB]
            ncompB = len(B.components)

            Pblock = P_density[A.col0:A.col0 + ncompA, B.col0:B.col0 + ncompB]
            if not np.any(Pblock):
                continue

            Ltot = A.L + B.L
            block = np.zeros((ncompA, ncompB, Nfield))

            for a in range(len(A.alphas)):
                for b in range(len(B.alphas)):
                    alpha = A.alphas[a]
                    beta = B.alphas[b]
                    p = alpha + beta
                    Pc = (alpha * A.center + beta * B.center) / p
                    cab = A.coeffs[a] * B.coeffs[b]
                    if cab == 0:
                        continue

                    Ex = _hermite_E_1d(alpha, beta, A.center[0], B.center[0], A.L, B.L)
                    Ey = _hermite_E_1d(alpha, beta, A.center[1], B.center[1], A.L, B.L)
                    Ez = _hermite_E_1d(alpha, beta, A.center[2], B.center[2], A.L, B.L)

                    PC = Pc - field_pts_bohr  # [Nfield, 3]
                    R = _hermite_R_full(p, PC, Ltot)  # [Ltot+1,Ltot+1,Ltot+1,Nfield]

                    prefac0 = (2 * np.pi / p) * cab

                    for ca in range(ncompA):
                        monoA = A.components[ca]["mono"]
                        coefA = A.components[ca]["coef"]
                        for cb in range(ncompB):
                            monoB = B.components[cb]["mono"]
                            coefB = B.components[cb]["coef"]
                            acc = np.zeros(Nfield)
                            for ta in range(monoA.shape[0]):
                                iA_, jA_, kA_ = monoA[ta]
                                Nprim_a = _norm_const_ijk(alpha, iA_, jA_, kA_)
                                for tb in range(monoB.shape[0]):
                                    iB_, jB_, kB_ = monoB[tb]
                                    w = coefA[ta] * coefB[tb]
                                    if w == 0:
                                        continue
                                    Nprim_b = _norm_const_ijk(beta, iB_, jB_, kB_)
                                    wterm = prefac0 * Nprim_a * Nprim_b * w
                                    for tt in range(0, iA_ + iB_ + 1):
                                        ex = Ex[iA_, iB_, tt]
                                        if ex == 0:
                                            continue
                                        for uu in range(0, jA_ + jB_ + 1):
                                            ey = Ey[jA_, jB_, uu]
                                            if ey == 0:
                                                continue
                                            for vv in range(0, kA_ + kB_ + 1):
                                                ez = Ez[kA_, kB_, vv]
                                                if ez == 0:
                                                    continue
                                                acc = acc + (wterm * ex * ey * ez) * R[tt, uu, vv, :]
                            block[ca, cb, :] += acc

            for ca in range(ncompA):
                for cb in range(ncompB):
                    if Pblock[ca, cb] == 0:
                        continue
                    V_elec += Pblock[ca, cb] * block[ca, cb, :]

    return V_elec
