function psi = g_eval_mo_on_grid(aobasis, mo_coeff, grid_pts_bohr)
%G_EVAL_MO_ON_GRID  Evaluates one or more MOs' real-space wavefunctions on
%   a grid of points [Ngrid x 3] (Bohr), accumulating shell-by-shell
%   directly into an [Ngrid x K] matrix -- never forming the full
%   [Ngrid x Nbasis] basis matrix, so memory stays O(Ngrid*K) regardless
%   of basis set size. Shared engine behind G_draw_mo_surface (K=1, a
%   single MO) and G_draw_density_surface (K=Nocc, all occupied MOs at
%   once, squared and occupation-weighted by the caller to build the
%   density -- evaluating every occupied MO in a single pass through the
%   shells this way is far cheaper than looping g_eval_mo_on_grid once
%   per occupied orbital, since the shell radial/angular parts, the
%   expensive part per grid point, are computed once and reused for
%   every column of MO_COEFF).
%
%   psi = g_eval_mo_on_grid(aobasis, mo_coeff, grid_pts_bohr)
%
%   mo_coeff: [Nbasis x K] MO coefficient columns (K=1 for a single MO).
%   Returns psi: [Ngrid x K].
%
%   Method: each contracted Gaussian-type basis function is evaluated
%   directly on the grid (primitive normalization N(alpha,l,m,n) =
%   (2*alpha/pi)^(3/4) * sqrt((4*alpha)^(l+m+n) / ((2l-1)!!(2m-1)!!(2n-1)!!));
%   pure D/F/G shells use the real solid harmonic polynomials of
%   Ribaldone & Desmarais (2024, arXiv:2412.16733, Table I -- a
%   re-derivation of Schlegel & Frisch, Int. J. Quantum Chem. 54, 83
%   (1995)), combined from individually-normalized Cartesian components
%   with the overall normalization computed at runtime (not hard-coded)
%   from the exact combinatorial self-overlap of same-center Cartesian
%   Gaussians -- see PURE_HARMONIC_VALUE/CART_COMPONENT_VALUE below.
%   Cartesian shell component ordering (including the reversed order used
%   natively by Gaussian .fchk files for L>=4) follows IOData's
%   (iodata.formats.fchk) documented Gaussian-native convention. Only S,
%   P, SP, D, F, G shells are handled -- call G_CHECK_SUPPORTED_SHELLS
%   first to raise a clear error on anything else, rather than silently
%   skipping unsupported shells here.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

    Ngrid  = size(grid_pts_bohr, 1);
    K      = size(mo_coeff, 2);
    Nshell = numel(aobasis.shell_types);
    psi = zeros(Ngrid, K);

    prim_offset = 0;
    col = 0;
    for s = 1:Nshell
        t     = aobasis.shell_types(s);
        nprim = aobasis.n_prim_per_shell(s);
        idx_p = prim_offset + (1:nprim);
        alphas = aobasis.prim_exponents(idx_p);
        cS     = aobasis.contraction_coeff(idx_p);
        center = aobasis.shell_coords_bohr(s, :);

        dx = grid_pts_bohr(:,1) - center(1);
        dy = grid_pts_bohr(:,2) - center(2);
        dz = grid_pts_bohr(:,3) - center(3);
        r2 = dx.^2 + dy.^2 + dz.^2;

        switch t
            case 0   % S
                radial = local_radial(alphas, cS, r2, 0);
                psi = psi + radial .* mo_coeff(col+1, :);
                col = col + 1;

            case 1   % P: X,Y,Z
                radial = local_radial(alphas, cS, r2, 1);
                psi = psi + radial .* (dx .* mo_coeff(col+1,:) + dy .* mo_coeff(col+2,:) + dz .* mo_coeff(col+3,:));
                col = col + 3;

            case -1  % SP: S, X, Y, Z
                cP = aobasis.sp_contraction_coeff(idx_p);
                radialS = local_radial(alphas, cS, r2, 0);
                radialP = local_radial(alphas, cP, r2, 1);
                psi = psi + radialS .* mo_coeff(col+1,:) ...
                           + radialP .* (dx .* mo_coeff(col+2,:) + dy .* mo_coeff(col+3,:) + dz .* mo_coeff(col+4,:));
                col = col + 4;

            case -2  % pure D: z2, xz, yz, x2-y2, xy
                radialD = local_radial(alphas, cS, r2, 2);
                Dxx = dx.*dx .* radialD;
                Dyy = dy.*dy .* radialD;
                Dzz = dz.*dz .* radialD;
                Dxy = dx.*dy .* radialD * sqrt(3);
                Dxz = dx.*dz .* radialD * sqrt(3);
                Dyz = dy.*dz .* radialD * sqrt(3);
                d_z2   = Dzz - 0.5*Dxx - 0.5*Dyy;
                d_x2y2 = (sqrt(3)/2) * (Dxx - Dyy);
                psi = psi + d_z2 .* mo_coeff(col+1,:) + Dxz .* mo_coeff(col+2,:) + Dyz .* mo_coeff(col+3,:) ...
                           + d_x2y2 .* mo_coeff(col+4,:) + Dxy .* mo_coeff(col+5,:);
                col = col + 5;

            case 2   % Cartesian D: xx,yy,zz,xy,xz,yz
                radialD = local_radial(alphas, cS, r2, 2);
                psi = psi + radialD .* ( ...
                    dx.*dx .* mo_coeff(col+1,:) + dy.*dy .* mo_coeff(col+2,:) + dz.*dz .* mo_coeff(col+3,:) ...
                  + sqrt(3)*(dx.*dy .* mo_coeff(col+4,:) + dx.*dz .* mo_coeff(col+5,:) + dy.*dz .* mo_coeff(col+6,:)) );
                col = col + 6;

            case -3  % pure F (7 components), Gaussian order m = 0,+1,-1,+2,-2,+3,-3
                radialF = local_radial(alphas, cS, r2, 3);
                terms = f_pure_terms();
                for c = 1:7
                    psi = psi + pure_harmonic_value(terms{c}, dx, dy, dz, radialF, 3) .* mo_coeff(col+c,:);
                end
                col = col + 7;

            case 3   % Cartesian F (10 components), Gaussian native order:
                     % xxx,yyy,zzz,xyy,xxy,xxz,xzz,yzz,yyz,xyz
                radialF = local_radial(alphas, cS, r2, 3);
                cart_f = [3 0 0; 0 3 0; 0 0 3; 1 2 0; 2 1 0; 2 0 1; 1 0 2; 0 1 2; 0 2 1; 1 1 1];
                for c = 1:10
                    tuv = cart_f(c,:);
                    psi = psi + cart_component_value(tuv(1), tuv(2), tuv(3), dx, dy, dz, radialF, 3) .* mo_coeff(col+c,:);
                end
                col = col + 10;

            case -4  % pure G (9 components), Gaussian order m = 0,+1,-1,+2,-2,+3,-3,+4,-4
                radialG = local_radial(alphas, cS, r2, 4);
                terms = g_pure_terms();
                for c = 1:9
                    psi = psi + pure_harmonic_value(terms{c}, dx, dy, dz, radialG, 4) .* mo_coeff(col+c,:);
                end
                col = col + 9;

            case 4   % Cartesian G (15 components), Gaussian native order (reversed
                     % lexicographic relative to the generic ascending-x convention):
                     % zzzz,yzzz,yyzz,yyyz,yyyy,xzzz,xyzz,xyyz,xyyy,xxzz,xxyz,xxyy,xxxz,xxxy,xxxx
                radialG = local_radial(alphas, cS, r2, 4);
                cart_g = [0 0 4; 0 1 3; 0 2 2; 0 3 1; 0 4 0; 1 0 3; 1 1 2; 1 2 1; 1 3 0; ...
                          2 0 2; 2 1 1; 2 2 0; 3 0 1; 3 1 0; 4 0 0];
                for c = 1:15
                    tuv = cart_g(c,:);
                    psi = psi + cart_component_value(tuv(1), tuv(2), tuv(3), dx, dy, dz, radialG, 4) .* mo_coeff(col+c,:);
                end
                col = col + 15;
        end

        prim_offset = prim_offset + nprim;
    end
end % g_eval_mo_on_grid


% =========================================================================
%  Local functions
% =========================================================================

function val = cart_component_value(t, u, v, dx, dy, dz, radial_shared, L)
%CART_COMPONENT_VALUE  Value of one Cartesian component (t,u,v), t+u+v=L,
%   of a Cartesian D/F/G shell, correcting RADIAL_SHARED (normalized for
%   the "diagonal", e.g. x^L, component) to this component's own
%   individual normalization: factor sqrt((2L-1)!! / ((2t-1)!!(2u-1)!!(2v-1)!!)).
    corr = sqrt(ddf(L) / (ddf(t)*ddf(u)*ddf(v)));
    val = corr * (dx.^t) .* (dy.^u) .* (dz.^v) .* radial_shared;
end

function val = pure_harmonic_value(terms, dx, dy, dz, radial_shared, L)
%PURE_HARMONIC_VALUE  Evaluates one normalized real solid harmonic pure
%   D/F/G component from its unnormalized polynomial TERMS = [P,t,u,v; ...]
%   (coefficients of x^t y^u z^v, from Ribaldone & Desmarais, a
%   re-derivation of Schlegel & Frisch, Int. J. Quantum Chem. 54, 83
%   (1995), Table I), combined with RADIAL_SHARED (normalized for the
%   "diagonal" x^L component). The overall normalization constant K is
%   computed here (not hard-coded) from the self-overlap of the
%   unnormalized combination, using the exact combinatorial formula for
%   the overlap of same-center, same-exponent Cartesian Gaussians:
%     <x^t1 y^u1 z^v1 | x^t2 y^u2 z^v2> (same alpha, normalized "diagonal"
%     radial part) = (t1+t2-1)!!(u1+u2-1)!!(v1+v2-1)!! / (2L-1)!!
    P = terms(:,1); t = terms(:,2); u = terms(:,3); v = terms(:,4);
    n = numel(P);
    S = 0;
    for k = 1:n
        for l = 1:n
            S = S + P(k)*P(l) * ddf((t(k)+t(l))/2) * ddf((u(k)+u(l))/2) * ddf((v(k)+v(l))/2);
        end
    end
    K = sqrt(ddf(L) / S);
    poly = zeros(size(dx));
    for k = 1:n
        poly = poly + P(k) * (dx.^t(k)) .* (dy.^u(k)) .* (dz.^v(k));
    end
    val = K * poly .* radial_shared;
end

function w = ddf(p)
%DDF  Odd double factorial (2p-1)!! indexed by nonnegative integer p,
%   with the standard convention ddf(0) = (-1)!! = 1.
    if p <= 0
        w = 1;
    else
        w = prod(1:2:(2*p-1));
    end
end

function terms = f_pure_terms()
%F_PURE_TERMS  Unnormalized real solid harmonic F polynomials (l=3),
%   Gaussian order m = 0,+1,-1,+2,-2,+3,-3. Each cell is [P,t,u,v; ...]
%   with P the coefficient of x^t y^u z^v. Source: Ribaldone & Desmarais
%   (2024), arXiv:2412.16733, Table I (re-derivation of Schlegel & Frisch
%   1995); cross-checked internally via the m=+3/-3 Re/Im[(x+iy)^3]
%   symmetry and, for l=2, against this function's own independently
%   derived D-shell formulas (exact match).
    terms = { ...
        [1 0 0 3; -1.5 2 0 1; -1.5 0 2 1]; ...   % m=0
        [6 1 0 2; -1.5 1 2 0; -1.5 3 0 0]; ...   % m=+1
        [6 0 1 2; -1.5 0 3 0; -1.5 2 1 0]; ...   % m=-1
        [-15 0 2 1; 15 2 0 1]; ...               % m=+2
        [30 1 1 1]; ...                          % m=-2
        [-45 1 2 0; 15 3 0 0]; ...                % m=+3
        [-15 0 3 0; 45 2 1 0] ...                 % m=-3
    };
end

function terms = g_pure_terms()
%G_PURE_TERMS  Unnormalized real solid harmonic G polynomials (l=4),
%   Gaussian order m = 0,+1,-1,+2,-2,+3,-3,+4,-4. See F_PURE_TERMS for
%   source/verification notes; m=+4/-4 cross-checked via the
%   Re/Im[(x+iy)^4] symmetry (both share the common prefactor 105).
    terms = { ...
        [1 0 0 4; -3 0 2 2; 0.375 0 4 0; -3 2 0 2; 0.75 2 2 0; 0.375 4 0 0]; ...  % m=0
        [10 1 0 3; -7.5 1 2 1; -7.5 3 0 1]; ...                                    % m=+1
        [10 0 1 3; -7.5 0 3 1; -7.5 2 1 1]; ...                                    % m=-1
        [-45 0 2 2; 7.5 0 4 0; 45 2 0 2; -7.5 4 0 0]; ...                          % m=+2
        [90 1 1 2; -15 1 3 0; -15 3 1 0]; ...                                      % m=-2
        [-315 1 2 1; 105 3 0 1]; ...                                               % m=+3
        [-105 0 3 1; 315 2 1 1]; ...                                               % m=-3
        [105 0 4 0; -630 2 2 0; 105 4 0 0]; ...                                    % m=+4
        [-420 1 3 0; 420 3 1 0] ...                                                % m=-4
    };
end

function radial = local_radial(alphas, coeffs, r2, L)
%LOCAL_RADIAL  Contracted radial part shared by all Cartesian components
%   of angular momentum L within one shell (using the "diagonal", e.g.
%   xx/x/S, primitive normalization):
%     N(alpha,L) = (2*alpha/pi)^(3/4) * sqrt((4*alpha)^L / (2L-1)!!)
%   (2L-1)!! for L=0..4 is 1,1,3,15,105 respectively.
    dfact = ddf(L);
    Nprim = (2*alphas/pi).^0.75 .* sqrt((4*alphas).^L / dfact);
    radial = zeros(size(r2));
    for pk = 1:numel(alphas)
        radial = radial + coeffs(pk) * Nprim(pk) * exp(-alphas(pk) * r2);
    end
end
