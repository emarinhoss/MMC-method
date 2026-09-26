function [Vl, El, VPl, EPl] = estimate_Vl(level_fn, L, N)
%ESTIMATE_VL  Convergence test: mean and variance of the MLMC corrections.
%
%   [Vl, El, VPl, EPl] = ESTIMATE_VL(level_fn, L, N)
%
%   Draws N samples on each level in the vector L (e.g. 0:5) with the
%   level sampler level_fn (same interface as for MLMC_ESTIMATE) and
%   returns, for every level, the mean and variance of the correction
%   Y_l = P_l - P_{l-1} and of P_l itself.  Each output has one row per
%   level and one column per output quantity.
%
%   Plotting log(Vl) and log(abs(El)) against l gives the strong and weak
%   convergence rates of the discretisation, and so shows whether MLMC is
%   worthwhile and which tolerance is realistic, before running
%   MLMC_ESTIMATE.  Example (tank problem, M = 4):
%
%       f = @(l,N) mlmc_tank_level(4, l, N, 1, 1, 4, 0, 1);
%       [Vl, El] = estimate_Vl(f, 0:5, 1e5);
%       semilogy(0:5, Vl, 'o-', 0:5, abs(El), 's-')
%
%   See also MLMC_ESTIMATE.

nl  = numel(L);
Vl  = [];  El  = [];  VPl = [];  EPl = [];
for i = 1:nl
    sums = level_fn(L(i), N);
    El(i,:)  = sums(1,:)/N;
    Vl(i,:)  = max(0, sums(2,:)/N - El(i,:).^2);
    EPl(i,:) = sums(3,:)/N;
    VPl(i,:) = max(0, sums(4,:)/N - EPl(i,:).^2);
end
end
