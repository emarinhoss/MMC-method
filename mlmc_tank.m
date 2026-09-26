%MLMC_TANK  MLMC for the tank-filling problem with a random outflow.
%
%   Estimates E[h(T)] for the height of liquid in a tank,
%       dh/dt = 10 + gamma*sin(t) - beta*sqrt(h),   h(0) = h0,
%       beta ~ Normal(beta_mean, beta_std^2),
%   using forward Euler on each level (see MLMC_TANK_LEVEL).  There is no
%   closed-form mean, so the result is compared with a reference value
%   obtained by Gauss-Hermite quadrature over beta, with each deterministic
%   ODE solved accurately by ode45.
%
%   See also MLMC_ESTIMATE, MLMC_TANK_LEVEL, TANK, DTANKFILL.

clear; close all; clc;
rng(0);                     % reproducible results

%% Problem parameters
T         = 1;              % final time
h0        = 1;              % initial height
gamma     = 4;              % amplitude of the periodic inflow
beta_mean = 0;              % mean of the outflow coefficient
beta_std  = 1;              % standard deviation of the outflow coefficient

%% MLMC parameters
M    = 4;                   % time-step refinement factor between levels
tol  = 1e-3;                % target RMS error
opts = struct('N0', 1e3, 'Lmin', 2, 'Lmax', 8, 'cost_exp', 1);

%% Multilevel Monte Carlo
disp('Multilevel Monte Carlo: tank filling with random outflow')
level_fn = @(l, N) mlmc_tank_level(M, l, N, T, h0, gamma, beta_mean, beta_std);
tic
[P, Nl, info] = mlmc_estimate(level_fn, M, tol, opts);
elapsed = toc;

%% Reference value: E[h(T)] by Gauss-Hermite quadrature over beta,
%  solving each deterministic ODE accurately with ode45.
%  Nodes xq and weights wq for the weight exp(-x^2) are computed with the
%  Golub-Welsch algorithm (eigen-decomposition of the Jacobi matrix).
nq = 20;
b  = sqrt((1:nq-1)/2);
[V, D] = eig(diag(b,1) + diag(b,-1));
[xq, iq] = sort(diag(D));
wq = sqrt(pi) * V(1,iq).'.^2;
ode_opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
href = 0;
for q = 1:nq
    beta = beta_mean + sqrt(2)*beta_std*xq(q);
    [~, h] = ode45(@(t,h) tankfill(t, h, gamma, beta), [0 T], h0, ode_opts);
    href = href + wq(q)*h(end);
end
href = href/sqrt(pi);

%% Report
fprintf('\n')
fprintf('  MLMC estimate      : %.6f\n', P)
fprintf('  reference (ode45)  : %.6f\n', href)
fprintf('  error              : %.2e   (tol = %g)\n', abs(P-href), tol)
fprintf('  finest level       : L = %d  (%d time steps)\n', info.L, M^info.L)
fprintf('  MLMC cost          : %.3g\n', info.cost)
fprintf('  standard MC cost   : %.3g  (estimated, same tol, finest level)\n', ...
        info.mc_cost)
fprintf('  savings factor     : %.1f\n', info.mc_cost/info.cost)
fprintf('  wall-clock time    : %.2f s\n', elapsed)

mlmc_plot(info, Nl, M, 'Tank filling with random outflow')

