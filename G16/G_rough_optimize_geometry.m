function xyzOut = G_rough_optimize_geometry(symbols, xyz, bondList, varargin)
%G_ROUGH_OPTIMIZE_GEOMETRY  Quick, approximate geometry clean-up for a
%   hand-built molecule.
%
%   xyzOut = G_ROUGH_OPTIMIZE_GEOMETRY(symbols, xyz, bondList)
%   xyzOut = G_ROUGH_OPTIMIZE_GEOMETRY(..., Name, Value, ...)
%
%   Pulls bond lengths toward the sum of covalent radii, spreads bond
%   angles toward a generic tetrahedral-like default, and pushes apart
%   non-bonded atoms that are clashing. This is explicitly NOT a real
%   force field (no MMFF/UFF atom typing, no torsional terms, no
%   electrostatics, no hybridization/aromaticity detection) and not a
%   substitute for a real Gaussian geometry optimization (the 'opt'
%   keyword) -- just enough to turn a rough, hand-placed sketch (e.g. from
%   G_MOLECULE_BUILDER) into a starting geometry that is not obviously
%   wrong (overlapping atoms, wildly off bond lengths) before sending it
%   to Gaussian for real.
%
%   Input:
%       symbols  - {Natoms x 1 cell} atomic symbols
%       xyz      - [Natoms x 3 double] starting Cartesian coordinates (Angstrom)
%       bondList - [Nbonds x 2] or [Nbonds x 3] Atom1/Atom2(/Order) pairs
%                  (same convention as G16_DRAW_MOLECULE's 'BondList')
%
%   Optional parameters (Name-Value):
%       'TargetAngleDeg' - generic target angle at every atom with 2+
%                          bonded neighbours, in degrees (default: 109.5,
%                          a generic tetrahedral-like default -- there is
%                          no hybridization/aromaticity detection here, so
%                          a genuinely planar/aromatic centre will not end
%                          up exactly right; refine those by hand
%                          afterward, e.g. by dragging in the 3D view)
%       'ClashFactor'    - a non-bonded pair is only penalised once its
%                          distance drops below ClashFactor times the sum
%                          of their covalent radii (default: 0.85 -- below
%                          the sum itself, so atoms that are merely close
%                          but not bonded are left alone)
%       'BondWeight', 'AngleWeight', 'ClashWeight' - relative weights
%                          between the three residual groups in the
%                          underlying least-squares fit (defaults: 1, 0.5, 2)
%
%   Minimises a sum-of-squares penalty over all atom coordinates via
%   LSQNONLIN (Optimization Toolbox). The penalty itself is invariant
%   under any rigid rotation/translation of the whole molecule, so without
%   correction the optimizer is free to also drift/spin it as a side
%   effect -- harmless for a Gaussian job, but visually disorienting right
%   after clicking "optimize"; the result is therefore re-centred on the
%   input's own centroid before being returned (rotation drift is left
%   alone -- a minor cosmetic point, not worth a full rigid-alignment step
%   for a "rough" clean-up).
%
%   A molecule with fewer than 2 atoms, or with no residual terms to
%   minimise at all (e.g. exactly 2 atoms with no bond declared between
%   them and no clash), is returned unchanged.
%
%   See also G_MOLECULE_BUILDER, G_ELEMENT_TABLE, LSQNONLIN.
%
%   Author: Sebastiano Trusso, CNR - Istituto per i Processi Chimico-Fisici (IPCF), Messina, Italy
%   Email: sebastiano.trusso@cnr.it
%   Developed with the assistance of an AI coding tool (Claude, Anthropic), under the author's supervision and review.

    p = inputParser;
    addRequired(p,  'symbols');
    addRequired(p,  'xyz',      @isnumeric);
    addRequired(p,  'bondList', @isnumeric);
    addParameter(p, 'TargetAngleDeg', 109.5, @isnumeric);
    addParameter(p, 'ClashFactor',    0.85,  @isnumeric);
    addParameter(p, 'BondWeight',     1.0,   @isnumeric);
    addParameter(p, 'AngleWeight',    0.5,   @isnumeric);
    addParameter(p, 'ClashWeight',    2.0,   @isnumeric);
    parse(p, symbols, xyz, bondList, varargin{:});
    opt = p.Results;

    symbols = cellstr(symbols(:));
    N = size(xyz, 1);
    if N < 2
        xyzOut = xyz;
        return
    end

    [~, covRadii] = G_element_table();
    radii = zeros(N, 1);
    for i = 1:N
        radii(i) = lookupRadius(covRadii, symbols{i});
    end

    bondedMask = false(N, N);
    nBonds = size(bondList, 1);
    targetLen = zeros(nBonds, 1);
    for b = 1:nBonds
        i = bondList(b,1); j = bondList(b,2);
        bondedMask(i,j) = true; bondedMask(j,i) = true;
        targetLen(b) = radii(i) + radii(j);
    end

    % Every pair of bonded neighbours sharing a common atom gets one angle
    % term, centred on that shared atom.
    angleTriples = zeros(0, 3);
    for i = 1:N
        nbrs = find(bondedMask(i,:));
        for a = 1:numel(nbrs)-1
            for c = a+1:numel(nbrs)
                angleTriples(end+1, :) = [i, nbrs(a), nbrs(c)]; %#ok<AGROW>
            end
        end
    end
    targetAngle = opt.TargetAngleDeg * pi / 180;

    nonBondedPairs = zeros(0, 3);
    for i = 1:N-1
        for j = i+1:N
            if ~bondedMask(i,j)
                nonBondedPairs(end+1, :) = [i, j, opt.ClashFactor * (radii(i) + radii(j))]; %#ok<AGROW>
            end
        end
    end

    if nBonds == 0 && isempty(angleTriples) && isempty(nonBondedPairs)
        xyzOut = xyz;
        return
    end

    x0 = xyz(:);
    % Levenberg-Marquardt explicitly: the default trust-region-reflective
    % algorithm requires at least as many residuals as free coordinates,
    % which a small/sparsely-connected molecule (e.g. very few non-bonded
    % pairs, few bonds) can easily have fewer of -- LSQNONLIN falls back
    % to Levenberg-Marquardt automatically in that case anyway (confirmed
    % during testing), so requesting it directly just avoids the warning.
    lsqOpts = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'levenberg-marquardt', ...
        'FunctionTolerance', 1e-8, 'StepTolerance', 1e-8, 'MaxFunctionEvaluations', 5000);
    try
        xSol = lsqnonlin(@residuals, x0, [], [], lsqOpts);
    catch ME
        error('G_rough_optimize_geometry:lsqnonlinFailed', 'Geometry clean-up failed: %s', ME.message);
    end
    xyzOut = reshape(xSol, N, 3);

    centroidOld = mean(xyz, 1);
    centroidNew = mean(xyzOut, 1);
    xyzOut = xyzOut + (centroidOld - centroidNew);

    % ---------------------------------------------------------------------
    function res = residuals(xflat)
        P = reshape(xflat, N, 3);

        rBond = zeros(nBonds, 1);
        for bi = 1:nBonds
            d = norm(P(bondList(bi,1),:) - P(bondList(bi,2),:));
            rBond(bi) = sqrt(opt.BondWeight) * (d - targetLen(bi));
        end

        nAngles = size(angleTriples, 1);
        rAngle = zeros(nAngles, 1);
        for ai = 1:nAngles
            ci = angleTriples(ai,1); n1 = angleTriples(ai,2); n2 = angleTriples(ai,3);
            v1 = P(n1,:) - P(ci,:);
            v2 = P(n2,:) - P(ci,:);
            cosTheta = dot(v1, v2) / (norm(v1) * norm(v2) + eps);
            cosTheta = max(min(cosTheta, 1), -1);
            theta = acos(cosTheta);
            rAngle(ai) = sqrt(opt.AngleWeight) * (theta - targetAngle);
        end

        nNB = size(nonBondedPairs, 1);
        rClash = zeros(nNB, 1);
        for qi = 1:nNB
            atomI = nonBondedPairs(qi,1); atomJ = nonBondedPairs(qi,2); thr = nonBondedPairs(qi,3);
            d = norm(P(atomI,:) - P(atomJ,:));
            if d < thr
                rClash(qi) = sqrt(opt.ClashWeight) * (thr - d);
            else
                rClash(qi) = 0;
            end
        end

        res = [rBond; rAngle; rClash];
    end
end

% -------------------------------------------------------------------------
function r = lookupRadius(covRadii, sym)
% Same fallback default (1.50 Angstrom) as G16_GET_BOND_LENGTH's own
% LOOKUP_RADIUS, for an element not in G_ELEMENT_TABLE's covalent-radii map.
    if isKey(covRadii, sym)
        r = covRadii(sym);
    else
        r = 1.50;
    end
end
