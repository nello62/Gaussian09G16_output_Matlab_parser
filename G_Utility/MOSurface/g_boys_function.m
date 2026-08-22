function F = g_boys_function(x, Nmax)
%G_BOYS_FUNCTION  Evaluates the Boys function F_n(x) = int_0^1 t^(2n) exp(-x t^2) dt
%   for n = 0..Nmax, for a vector of arguments x >= 0, via the exact closed
%   form in terms of the (regularized) lower incomplete gamma function
%       F_n(x) = gammainc(x, n+1/2) * gamma(n+1/2) / (2 x^(n+1/2)),  x > 0
%       F_n(0) = 1/(2n+1)
%   (obtained by the substitution u = x t^2 in the defining integral).
%   MATLAB's built-in GAMMAINC is a well-tested, numerically robust special
%   function valid over the full x range, so this avoids implementing a
%   separate small-x series / large-x asymptotic / recursion-direction
%   strategy by hand. Verified against direct numerical quadrature of the
%   defining integral to machine precision (relative error ~1e-15) for
%   x in [0.001, 60] and n in [0, 4].
%
%   F = g_boys_function(x, Nmax)
%
%   x: [Npts x 1] (or any shape), Nmax: scalar. Returns F: [Npts x (Nmax+1)],
%   column n+1 = F_n(x).
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

x = x(:);
Npts = numel(x);
F = zeros(Npts, Nmax+1);
small = x < 1e-12;

for n = 0:Nmax
    col = zeros(Npts,1);
    if any(~small)
        xv = x(~small);
        col(~small) = gammainc(xv, n+0.5) .* gamma(n+0.5) ./ (2*xv.^(n+0.5));
    end
    if any(small)
        col(small) = 1/(2*n+1);
    end
    F(:,n+1) = col;
end

end % g_boys_function
