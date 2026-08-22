function V_elec = g_eval_esp_analytic(aobasis, P_density, field_pts_bohr, verbose)
%G_EVAL_ESP_ANALYTIC  Exact analytic electronic-Coulomb-potential evaluation
%   V_elec(r) = sum_munu P_munu <phi_mu | 1/|r-r'| | phi_nu>
%   at each field point r, via the McMurchie-Davidson Hermite-Gaussian
%   scheme (E_t^{ij} expansion coefficients + R_tuv Hermite Coulomb
%   integrals + the Boys function), rather than the numerical grid-sum
%   approximation used by G_DRAW_ESP_SURFACE. See Theory_Vibrations_
%   Polarizability.tex, Part III, for the full derivation and the formula
%   sources (Helgaker, McMurchie & Davidson 1978).
%
%   V_elec = g_eval_esp_analytic(aobasis, P_density, field_pts_bohr, verbose)
%
%   aobasis: as read by G_READ_AOBASIS_FROM_FCHK. P_density: [Nbasis x
%   Nbasis] density matrix (e.g. 2*occ_coeff*occ_coeff' for a closed-shell
%   calculation). field_pts_bohr: [Nfield x 3].
%
%   ONLY S, P, SP, and pure D shells are supported (matching the basis
%   actually used for validation, see below); Cartesian D and F/G shells
%   raise a clear error rather than a silently wrong result. This is a
%   narrower scope than G_EVAL_MO_ON_GRID (S..G): analytic integrals need
%   substantially more machinery per angular-momentum level than simple
%   point evaluation, so this was scoped to the common case first.
%
%   Validated (not just derived) before use: the core E_t^{ij}/R_tuv
%   recursion machinery was checked against independent references at
%   each step -- an S-S overlap and a D-D 1D overlap from the E-recursion
%   against their closed-form/numerically-integrated values (relative
%   error ~1e-16); an S-S, a P-S, and a Px-Py (a genuine multi-index t,u
%   case) nuclear-attraction integral from the full E+R combination
%   against direct 3D numerical quadrature (relative error ~1e-9, the
%   quadrature's own precision limit); the vectorized (multi-field-point)
%   R_tuv builder against the single-point recursive version, exact match
%   (0.0 difference) for every (t,u,v) up to a D-D pair's combined degree.
%   The complete function is further validated end-to-end against a real
%   Gaussian-generated ESP cube file -- see the "Accuracy" note in
%   G_DRAW_ESP_SURFACE and Theory_Vibrations_Polarizability.tex Part III.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

if nargin < 4
    verbose = false;
end

Nfield = size(field_pts_bohr, 1);
V_elec = zeros(Nfield, 1);
Nshell = numel(aobasis.shell_types);

% -------------------------------------------------------------------------
% Expand each shell into one or more "subshells" (an SP shell becomes an
% independent S part, using .contraction_coeff, and P part, using
% .sp_contraction_coeff -- exactly as G_EVAL_MO_ON_GRID treats SP), each
% with a fixed angular momentum L and a list of final basis-function
% components, each a linear combination (already correctly normalized,
% via individually-normalized Cartesian monomials -- see the K,P_k
% derivation in Theory_Vibrations_Polarizability.tex Part III) of
% Cartesian monomials (l,m,n), l+m+n=L.
% -------------------------------------------------------------------------
subshells = struct('L', {}, 'alphas', {}, 'coeffs', {}, 'center', {}, 'components', {}, 'col0', {});
prim_offset = 0;
col = 0;
for s = 1:Nshell
    t = aobasis.shell_types(s);
    nprim = aobasis.n_prim_per_shell(s);
    idx_p = prim_offset + (1:nprim);
    alphas = aobasis.prim_exponents(idx_p);
    cS = aobasis.contraction_coeff(idx_p);
    center = aobasis.shell_coords_bohr(s, :);

    switch t
        case 0
            subshells(end+1) = make_subshell(0, alphas, cS, center, s_components(), col); %#ok<AGROW>
            col = col + 1;
        case 1
            subshells(end+1) = make_subshell(1, alphas, cS, center, p_components(), col); %#ok<AGROW>
            col = col + 3;
        case -1
            cP = aobasis.sp_contraction_coeff(idx_p);
            subshells(end+1) = make_subshell(0, alphas, cS, center, s_components(), col); %#ok<AGROW>
            col = col + 1;
            subshells(end+1) = make_subshell(1, alphas, cP, center, p_components(), col); %#ok<AGROW>
            col = col + 3;
        case -2
            subshells(end+1) = make_subshell(2, alphas, cS, center, pure_d_components(), col); %#ok<AGROW>
            col = col + 5;
        otherwise
            error('g_eval_esp_analytic:unsupportedShell', ...
                ['g_eval_esp_analytic only supports S, P, SP, and pure D shells ' ...
                 '(shell type %d at shell %d is not one of these). Cartesian D and ' ...
                 'F/G shells are not yet implemented for the analytic ESP evaluator ' ...
                 '(unlike g_eval_mo_on_grid, which supports up to G) -- analytic ' ...
                 'integrals need substantially more machinery per angular-momentum ' ...
                 'level than point evaluation.'], t, s);
    end
    prim_offset = prim_offset + nprim;
end

if verbose
    fprintf('g_eval_esp_analytic: %d subshells, %d basis functions, %d field points\n', ...
        numel(subshells), col, Nfield);
end

% -------------------------------------------------------------------------
% Shell-pair loop
% -------------------------------------------------------------------------
Nsub = numel(subshells);
tic;
for iA = 1:Nsub
    A = subshells(iA);
    ncompA = numel(A.components);
    for iB = 1:Nsub
        B = subshells(iB);
        ncompB = numel(B.components);

        Pblock = P_density(A.col0+(1:ncompA), B.col0+(1:ncompB));
        if all(Pblock(:) == 0)
            continue
        end

        Ltot = A.L + B.L;
        block = zeros(ncompA, ncompB, Nfield);

        for a = 1:numel(A.alphas)
            for b = 1:numel(B.alphas)
                alpha = A.alphas(a); beta = B.alphas(b);
                p = alpha + beta;
                Pc = (alpha*A.center + beta*B.center) / p;
                cab = A.coeffs(a) * B.coeffs(b);
                if cab == 0, continue; end

                % NOTE: the Gaussian-product prefactor exp(-mu*|A-B|^2) is
                % already carried through the ENTIRE hermite_E_1d recursion
                % (it is baked into the E(1,1,1) seed and propagates
                % linearly to every (t,u,v) entry), since the product of
                % the three 1D seeds already equals the full 3D prefactor.
                % Do NOT multiply by it again here (that was Bug #3: it
                % squared the pair-distance suppression for every
                % different-center primitive pair, invisible for
                % same-center terms where the factor is 1).
                Ex = hermite_E_1d(alpha, beta, A.center(1), B.center(1), A.L, B.L);
                Ey = hermite_E_1d(alpha, beta, A.center(2), B.center(2), A.L, B.L);
                Ez = hermite_E_1d(alpha, beta, A.center(3), B.center(3), A.L, B.L);

                PC = Pc - field_pts_bohr;   % [Nfield x 3]
                R = hermite_R_full(p, PC, Ltot);   % [Ltot+1,Ltot+1,Ltot+1,Nfield]

                prefac0 = (2*pi/p) * cab;

                for ca = 1:ncompA
                    monoA = A.components(ca).mono; coefA = A.components(ca).coef;
                    for cb = 1:ncompB
                        monoB = B.components(cb).mono; coefB = B.components(cb).coef;
                        acc = zeros(Nfield,1);
                        for ta = 1:size(monoA,1)
                            iA_ = monoA(ta,1); jA_ = monoA(ta,2); kA_ = monoA(ta,3);
                            % Individual (not shared-diagonal) primitive
                            % normalization for THIS specific monomial --
                            % off-diagonal monomials (e.g. xy, xz, yz)
                            % need a different constant from diagonal ones
                            % (e.g. xx) at the same total degree L; see
                            % the header note on the bug this fixes.
                            Nprim_a = norm_const_ijk(alpha, iA_, jA_, kA_);
                            for tb = 1:size(monoB,1)
                                iB_ = monoB(tb,1); jB_ = monoB(tb,2); kB_ = monoB(tb,3);
                                w = coefA(ta)*coefB(tb);
                                if w == 0, continue; end
                                Nprim_b = norm_const_ijk(beta, iB_, jB_, kB_);
                                wterm = prefac0 * Nprim_a * Nprim_b * w;
                                for tt = 0:(iA_+iB_)
                                    ex = Ex(iA_+1,iB_+1,tt+1);
                                    if ex == 0, continue; end
                                    for uu = 0:(jA_+jB_)
                                        ey = Ey(jA_+1,jB_+1,uu+1);
                                        if ey == 0, continue; end
                                        for vv = 0:(kA_+kB_)
                                            ez = Ez(kA_+1,kB_+1,vv+1);
                                            if ez == 0, continue; end
                                            acc = acc + (wterm*ex*ey*ez) * reshape(R(tt+1,uu+1,vv+1,:),[],1);
                                        end
                                    end
                                end
                            end
                        end
                        block(ca,cb,:) = reshape(block(ca,cb,:),[],1) + acc;
                    end
                end
            end
        end

        for ca = 1:ncompA
            for cb = 1:ncompB
                if Pblock(ca,cb) == 0, continue; end
                V_elec = V_elec + Pblock(ca,cb) * reshape(block(ca,cb,:),[],1);
            end
        end
    end
    if verbose && mod(iA, max(1,round(Nsub/10))) == 0
        fprintf('  shell-pair loop: %d/%d subshells (%.1f s elapsed)\n', iA, Nsub, toc);
    end
end

end % g_eval_esp_analytic


% =========================================================================
%  Local functions
% =========================================================================

function sub = make_subshell(L, alphas, coeffs, center, components, col0)
    sub.L = L;
    sub.alphas = alphas;
    sub.coeffs = coeffs;
    sub.center = center;
    sub.components = components;
    sub.col0 = col0;
end

function c = s_components()
    c = struct('mono', {[0 0 0]}, 'coef', {1});
end

function c = p_components()
    c(1) = struct('mono', [1 0 0], 'coef', 1);
    c(2) = struct('mono', [0 1 0], 'coef', 1);
    c(3) = struct('mono', [0 0 1], 'coef', 1);
end

function c = pure_d_components()
%PURE_D_COMPONENTS  The 5 pure-D basis functions (Gaussian order m =
%   0,+1,-1,+2,-2), each as a linear combination of individually-
%   normalized Cartesian D monomials (order: xx,yy,zz,xy,xz,yz), with the
%   overall normalization K folded into the returned coefficients.
%
%   K uses the TRUE per-term-pair individually-normalized overlap
%       <mono_k|mono_l> = raw(k,l) / sqrt(D(k)*D(l)),   D(m) = ddf(m1)ddf(m2)ddf(m3)
%   NOT raw(k,l)/(2L-1)!! (a single shared divisor) -- an earlier version
%   used the latter, which happens to coincide with the former whenever
%   every term in a component has D(m)=(2L-1)!! (true for z2 and x2-y2,
%   whose terms are all "diagonal" xx/yy/zz-type), but is WRONG whenever
%   D(m) differs from (2L-1)!! -- true for any single-term off-diagonal
%   component (xz, yz, xy, each with D=1 vs (2L-1)!!=3), which the buggy
%   version overnormalized by a factor of sqrt(3), giving a self-overlap
%   of 3 instead of 1. Caught by an analytic-overlap self-consistency
%   check (diag(S) should be all-ones for properly normalized basis
%   functions; it was not) after an end-to-end ESP comparison against a
%   real Gaussian cube file first revealed a several-percent deficit in
%   the analytic electronic potential -- see G_EVAL_ESP_ANALYTIC's header.
    terms = { ...
        [1 0 0 2; -0.5 2 0 0; -0.5 0 2 0]; ...   % m=0 (z2): zz - 0.5xx - 0.5yy
        [1 1 0 1]; ...                           % m=+1 (xz)
        [1 0 1 1]; ...                           % m=-1 (yz)
        [0.5*sqrt(3) 2 0 0; -0.5*sqrt(3) 0 2 0]; ... % m=+2 (x2-y2): (sqrt3/2)(xx-yy)
        [1 1 1 0] ...                            % m=-2 (xy)
    };
    c = repmat(struct('mono',[],'coef',[]), 1, 5);
    for k = 1:5
        Tk = terms{k};
        P = Tk(:,1); mono = Tk(:,2:4);
        nterm = size(Tk,1);
        D = zeros(nterm,1);
        for i = 1:nterm
            D(i) = ddf(mono(i,1))*ddf(mono(i,2))*ddf(mono(i,3));
        end
        S = 0;
        for i = 1:nterm
            for j = 1:nterm
                raw = ddf_pair(mono(i,1),mono(j,1)) * ddf_pair(mono(i,2),mono(j,2)) * ddf_pair(mono(i,3),mono(j,3));
                S = S + P(i)*P(j) * raw / sqrt(D(i)*D(j));
            end
        end
        K = 1/sqrt(S);
        c(k).mono = mono;
        c(k).coef = K * P;
    end
end

function v = ddf_pair(a,b)
    v = ddf((a+b)/2);
end

function w = ddf(p)
    if p <= 0
        w = 1;
    else
        w = prod(1:2:(2*p-1));
    end
end

function N = norm_const_ijk(alpha, i, j, k)
%NORM_CONST_IJK  TRUE individual primitive normalization for the specific
%   Cartesian monomial x^i y^j z^k, Eq. (gto-norm) of
%   Theory_Vibrations_Polarizability.tex Part III:
%     N(alpha,i,j,k) = (2alpha/pi)^(3/4) sqrt((4alpha)^L / (ddf(i)ddf(j)ddf(k)))
%   Must be evaluated per MONOMIAL, not once per shell from L alone: an
%   earlier version used a single shared N(alpha,L) (the "diagonal", e.g.
%   x^L, value) for every monomial in a pure-D component -- correct for
%   the diagonal monomials (xx,yy,zz) but silently wrong (undernormalized
%   by a factor of sqrt(3)) for the off-diagonal ones (xy,xz,yz), since
%   ddf(1)*ddf(1)*ddf(0)=1 there rather than ddf(2)=3. Caught by an
%   end-to-end comparison against a real Gaussian-generated ESP cube file
%   (4-NTP): analytic V_elec was ~6-7% low at a point dominated by a
%   pure-D shell pair, which is a large error specifically because V_nuc
%   and V_elec nearly cancel outside the molecule (see G_DRAW_ESP_SURFACE
%   and Theory_Vibrations_Polarizability.tex Part III for the physical
%   reason this cancellation makes even a small V_elec error dangerous).
    L = i+j+k;
    N = (2*alpha/pi)^0.75 * sqrt((4*alpha)^L / (ddf(i)*ddf(j)*ddf(k)));
end

function E = hermite_E_1d(alpha, beta, AX, BX, La, Lb)
%HERMITE_E_1D  McMurchie-Davidson Hermite expansion coefficients E_t^{ij}
%   for one Cartesian direction, i=0..La, j=0..Lb, t=0..i+j. Verified
%   against independent overlap-integral references, see this file's
%   header and G_EVAL_ESP_ANALYTIC's help text.
    p = alpha+beta; mu = alpha*beta/p; P = (alpha*AX+beta*BX)/p;
    XPA = P-AX; XPB = P-BX; XAB = AX-BX;
    Tmax = La+Lb;
    E = zeros(La+1, Lb+1, Tmax+1);
    E(1,1,1) = exp(-mu*XAB^2);
    for j = 1:Lb
        ij_prev = j-1;
        for t = 0:j
            term1 = 0; if t>=1 && (t-1)<=ij_prev, term1 = E(1,j,t); end
            term2 = 0; if t<=ij_prev, term2 = E(1,j,t+1); end
            term3 = 0; if (t+1)<=ij_prev, term3 = E(1,j,t+2); end
            E(1,j+1,t+1) = term1/(2*p) + XPB*term2 + (t+1)*term3;
        end
    end
    for i = 1:La
        for j = 0:Lb
            ij_prev = (i-1)+j;
            for t = 0:(i+j)
                term1 = 0; if t>=1 && (t-1)<=ij_prev, term1 = E(i,j+1,t); end
                term2 = 0; if t<=ij_prev, term2 = E(i,j+1,t+1); end
                term3 = 0; if (t+1)<=ij_prev, term3 = E(i,j+1,t+2); end
                E(i+1,j+1,t+1) = term1/(2*p) + XPA*term2 + (t+1)*term3;
            end
        end
    end
end

function R = hermite_R_full(p, PC, Ltot)
%HERMITE_R_FULL  McMurchie-Davidson Hermite Coulomb integrals R_tuv (n=0),
%   for all t+u+v<=Ltot, vectorized over field points. PC: [Nfield x 3]
%   (P - C for each field point C). Verified to exactly match (0.0
%   difference) an independent single-point recursive implementation for
%   every (t,u,v) up to a D-D pair's combined degree -- see this file's
%   header and G_EVAL_ESP_ANALYTIC's help text.
    Nfield = size(PC,1);
    Rfull = zeros(Ltot+1,Ltot+1,Ltot+1,Ltot+1,Nfield);

    x = p*sum(PC.^2,2);
    Fn = g_boys_function(x, Ltot);
    for n = 0:Ltot
        Rfull(1,1,1,n+1,:) = reshape(((-2*p)^n) * Fn(:,n+1), 1,1,1,1,Nfield);
    end

    for D = 1:Ltot
        for t = 0:D
            for u = 0:(D-t)
                v = D-t-u;
                if t > 0
                    for n = 0:(Ltot-D)
                        term = PC(:,1) .* reshape(Rfull(t,u+1,v+1,n+2,:),[],1);
                        if t >= 2
                            term = term + (t-1)*reshape(Rfull(t-1,u+1,v+1,n+2,:),[],1);
                        end
                        Rfull(t+1,u+1,v+1,n+1,:) = reshape(term,1,1,1,1,Nfield);
                    end
                elseif u > 0
                    for n = 0:(Ltot-D)
                        term = PC(:,2) .* reshape(Rfull(1,u,v+1,n+2,:),[],1);
                        if u >= 2
                            term = term + (u-1)*reshape(Rfull(1,u-1,v+1,n+2,:),[],1);
                        end
                        Rfull(1,u+1,v+1,n+1,:) = reshape(term,1,1,1,1,Nfield);
                    end
                else
                    for n = 0:(Ltot-D)
                        term = PC(:,3) .* reshape(Rfull(1,1,v,n+2,:),[],1);
                        if v >= 2
                            term = term + (v-1)*reshape(Rfull(1,1,v-1,n+2,:),[],1);
                        end
                        Rfull(1,1,v+1,n+1,:) = reshape(term,1,1,1,1,Nfield);
                    end
                end
            end
        end
    end

    R = reshape(Rfull(:,:,:,1,:), Ltot+1, Ltot+1, Ltot+1, Nfield);
end
