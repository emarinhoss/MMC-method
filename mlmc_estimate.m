function [P, Nl, info] = mlmc_estimate(level_fn, M, tol, opts)
%MLMC_ESTIMATE  Adaptive Multilevel Monte Carlo (MLMC) estimator.
%
%   [P, Nl, info] = MLMC_ESTIMATE(level_fn, M, tol)
%   [P, Nl, info] = MLMC_ESTIMATE(level_fn, M, tol, opts)
%
%   Estimates E[P] for a quantity of interest P that is computed with a
%   discretisation (time step and/or mesh size) refined by a factor M from
%   one level to the next.  The algorithm follows
%
%     M.B. Giles, "Multilevel Monte Carlo path simulation",
%     Operations Research 56(3):607-617, 2008.
%
%   and uses the telescoping sum
%
%     E[P_L] = E[P_0] + sum_{l=1}^{L} E[P_l - P_{l-1}],
%
%   where every correction Y_l = P_l - P_{l-1} is estimated independently
%   with N_l samples.  Levels are added until the estimated bias is below
%   tol/sqrt(2), and the N_l are chosen to make the sampling variance
%   tol^2/2 at minimum cost, so that the total root-mean-square error is
%   approximately tol.
%
%   Inputs
%     level_fn  function handle, sums = level_fn(l, N).  For level l it must
%               draw N independent samples and return a 4-by-Q array
%                 sums(1,:) = sum of Y_l          (Y_0 = P_0)
%                 sums(2,:) = sum of Y_l.^2
%                 sums(3,:) = sum of P_l
%                 sums(4,:) = sum of P_l.^2
%               Q is the number of output quantities (Q = 1 for a scalar,
%               Q > 1 e.g. for a solution sampled at several grid points).
%               For vector outputs the worst case over the Q components is
%               used for the variance and bias estimates.
%     M         refinement factor between consecutive levels (M >= 2).
%     tol       target root-mean-square error of the estimate.
%     opts      optional struct with fields (defaults in brackets)
%                 N0       initial number of samples on a new level [1e4]
%                 Lmin     minimum level before testing convergence  [2]
%                 Lmax     maximum level; stop with a warning        [10]
%                 cost_exp cost of one sample on level l is taken as
%                          M^(cost_exp*l) (1 for an ODE, 2 for a 1D PDE
%                          refined in both space and time)            [1]
%                 verbose  print progress                           [true]
%
%   Outputs
%     P     MLMC estimate of E[P] (1-by-Q).
%     Nl    number of samples actually used on each level (1-by-(L+1)).
%     info  struct with per-level diagnostics:
%             L      finest level used
%             El     mean of Y_l          ((L+1)-by-Q)
%             Vl     variance of Y_l      ((L+1)-by-Q)
%             EPl    mean of P_l          ((L+1)-by-Q)
%             VPl    variance of P_l      ((L+1)-by-Q)
%             Cl     cost per sample on each level (1-by-(L+1))
%             cost   total MLMC cost, sum(Nl.*Cl)
%             mc_cost  estimated cost of standard Monte Carlo on the finest
%                      level for the same tolerance, 2*max(VPl(end,:))/tol^2*Cl(end)
%             converged  false if Lmax was reached before convergence
%
%   See also ESTIMATE_VL, MLMC_PLOT.

if nargin < 4, opts = struct(); end
N0       = get_opt(opts, 'N0',       1e4);
Lmin     = get_opt(opts, 'Lmin',     2);
Lmax     = get_opt(opts, 'Lmax',     10);
cost_exp = get_opt(opts, 'cost_exp', 1);
verbose  = get_opt(opts, 'verbose',  true);

L     = -1;
suml  = [];      % suml(:,:,l+1) accumulates the 4-by-Q sums on level l
Nl    = [];      % samples taken so far on each level
converged = false;

while ~converged
    %
    % Step 1: add a new level and take an initial set of samples on it,
    %         so that its variance can be estimated.
    %
    L = L + 1;
    sums = level_fn(L, N0);
    suml(:,:,L+1) = sums;
    Nl(L+1) = N0;

    %
    % Step 2: estimate the variance V_l of Y_l on every level.
    %         V(Y) = E[Y^2] - E[Y]^2, clipped at 0 to avoid round-off.
    %
    [El, Vl] = level_moments(suml, Nl, 1);
    V  = max(Vl, [], 2).';               % worst case over the Q outputs
    Cl = M.^(cost_exp*(0:L));            % cost of one sample per level

    %
    % Step 3: optimal number of samples per level (Giles 2008):
    %           N_l = ceil( 2/tol^2 * sqrt(V_l/C_l) * sum_k sqrt(V_k*C_k) )
    %         which minimises the total cost subject to a sampling
    %         variance of tol^2/2.
    %
    Nopt = ceil(2/tol^2 * sqrt(V./Cl) * sum(sqrt(V.*Cl)));

    %
    % Step 4: draw the extra samples required on each level.
    %
    for l = 0:L
        dNl = Nopt(l+1) - Nl(l+1);
        if dNl > 0
            suml(:,:,l+1) = suml(:,:,l+1) + level_fn(l, dNl);
            Nl(l+1) = Nl(l+1) + dNl;
        end
    end

    if verbose
        fprintf('  L = %2d   Nl = %s\n', L, sprintf('%d ', Nl));
    end

    %
    % Step 5: convergence test.  Assuming first order weak convergence,
    %         E[Y_l] ~ M^(-l), so the remaining bias is ~ E[Y_L]/(M-1).
    %         Both of the last two corrections are used (the older one
    %         scaled by 1/M) to make the test more robust (Giles 2008).
    %
    if L >= Lmin
        El = level_moments(suml, Nl, 1);
        bias = max(max(abs(El(L,:))/M), max(abs(El(L+1,:))));
        converged = bias < (M-1)*tol/sqrt(2);
    end

    if ~converged && L >= Lmax
        warning('mlmc_estimate:Lmax', ...
            'Reached Lmax = %d before the bias test was satisfied.', Lmax);
        break
    end
end

%
% Evaluate the multilevel estimator  P = sum_l mean(Y_l).
%
[El, Vl]   = level_moments(suml, Nl, 1);
[EPl, VPl] = level_moments(suml, Nl, 3);
P = sum(El, 1);

Cl = M.^(cost_exp*(0:L));
info = struct('L', L, 'El', El, 'Vl', Vl, 'EPl', EPl, 'VPl', VPl, ...
              'Cl', Cl, 'cost', sum(Nl.*Cl), ...
              'mc_cost', 2*max(VPl(end,:))/tol^2*Cl(end), ...
              'converged', converged);
end

% -------------------------------------------------------------------------
function [E, V] = level_moments(suml, Nl, row)
%LEVEL_MOMENTS  Sample mean and variance on each level from the running
%   sums stored in rows ROW (sum) and ROW+1 (sum of squares) of SUML.
nl = numel(Nl);
Q  = size(suml, 2);
E  = zeros(nl, Q);
V  = zeros(nl, Q);
for l = 1:nl
    E(l,:) = suml(row,  :,l) / Nl(l);
    V(l,:) = max(0, suml(row+1,:,l) / Nl(l) - E(l,:).^2);
end
end

function v = get_opt(opts, name, default)
%GET_OPT  Field NAME of OPTS, or DEFAULT if it is absent.
if isfield(opts, name)
    v = opts.(name);
else
    v = default;
end
end
