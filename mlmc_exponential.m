%MLMC_EXPONENTIAL  MLMC for exponential growth with a random rate.
%
%   Estimates E[y(T)] for
%       dy/dt = k*y,   y(0) = y0,   k ~ Normal(k_mean, k_std^2),
%   using forward Euler on each level (see MLMC_EXPONENTIAL_LEVEL) and
%   compares the result with the exact mean
%       E[y(T)] = y0 * exp(k_mean*T + (k_std*T)^2/2),
%   which follows from the moment generating function of a normal random
%   variable.  (Note that this is larger than y0*exp(k_mean*T), the
%   solution for the mean rate.)
%
%   See also MLMC_ESTIMATE, MLMC_EXPONENTIAL_LEVEL, MLMC_PLOT.

clear; close all; clc;
rng(0);                     % reproducible results

%% Problem parameters
y0     = 10;                % initial condition
T      = 1;                 % final time
k_mean = 3;                 % mean growth rate
k_std  = 0.5;               % standard deviation of the growth rate

%% MLMC parameters
M    = 3;                   % time-step refinement factor between levels
tol  = 0.5;                 % target RMS error (the answer is ~228)
opts = struct('N0', 1e3, ...       % initial samples on each new level
              'Lmin', 3, ...       % do not test convergence before this
              'Lmax', 10, ...      % give up after this level
              'cost_exp', 1);      % cost per sample ~ number of steps

%% Exact solution
exact = y0*exp(k_mean*T + 0.5*(k_std*T)^2);

%% Multilevel Monte Carlo
disp('Multilevel Monte Carlo: exponential growth with random rate')
level_fn = @(l, N) mlmc_exponential_level(M, l, N, T, y0, k_mean, k_std);
tic
[P, Nl, info] = mlmc_estimate(level_fn, M, tol, opts);
elapsed = toc;

%% Report
fprintf('\n')
fprintf('  MLMC estimate      : %.4f\n', P)
fprintf('  exact mean         : %.4f\n', exact)
fprintf('  error              : %.4f   (tol = %g)\n', abs(P-exact), tol)
fprintf('  finest level       : L = %d  (%d time steps)\n', info.L, M^info.L)
fprintf('  MLMC cost          : %.3g\n', info.cost)
fprintf('  standard MC cost   : %.3g  (estimated, same tol, finest level)\n', ...
        info.mc_cost)
fprintf('  savings factor     : %.1f\n', info.mc_cost/info.cost)
fprintf('  wall-clock time    : %.2f s\n', elapsed)

mlmc_plot(info, Nl, M, 'Exponential growth with random rate')
