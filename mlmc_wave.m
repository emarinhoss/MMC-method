%MLMC_WAVE  MLMC for the 1D advection (one-way wave) equation with a random
%   wave speed.
%
%   Estimates the mean solution E[u(x,T)] of
%       u_t + c*u_x = 0,   x in [-R, R],   u(x,0) = exp(-x^2),
%       c ~ Normal(c_mean, c_std^2),
%   at the grid points of the coarsest mesh, using first-order upwind
%   finite differences in which the mesh size and the time step are both
%   halved from one level to the next (see MLMC_WAVE_LEVEL).
%
%   Since u(x,T) = exp(-(x - c*T)^2) and c*T ~ Normal(mu, s^2) with
%   mu = c_mean*T and s = c_std*T, the exact mean is
%       E[u(x,T)] = exp(-(x - mu)^2 / (1 + 2*s^2)) / sqrt(1 + 2*s^2).
%
%   See also MLMC_ESTIMATE, MLMC_WAVE_LEVEL, MLMC_PLOT.

clear; close all; clc;
rng(0);                     % reproducible results

%% Problem parameters
T      = 1;                 % final time
R      = 5;                 % half-width of the domain [-R, R]
nr     = 64;                % number of cells on the coarsest grid (level 0)
c_mean = -1;                % mean wave speed
c_std  = 0.2;               % standard deviation of the wave speed
cfl    = 0.9;               % CFL number used to choose the time step

%% MLMC parameters
M    = 2;                   % refinement factor (fixed by MLMC_WAVE_LEVEL)
tol  = 2e-3;                % target RMS error (max over the grid points)
opts = struct('N0', 100, ...       % initial samples on each new level
              'Lmin', 2, ...
              'Lmax', 8, ...
              'cost_exp', 2);      % space AND time refined: cost ~ 4^l

%% Multilevel Monte Carlo
disp('Multilevel Monte Carlo: 1D advection with random wave speed')
level_fn = @(l, N) mlmc_wave_level(l, N, T, R, nr, c_mean, c_std, cfl);
tic
[P, Nl, info] = mlmc_estimate(level_fn, M, tol, opts);
elapsed = toc;

%% Exact mean solution
x     = linspace(-R, R, nr+1);
mu    = c_mean*T;
s2    = (c_std*T)^2;
exact = exp(-(x - mu).^2 / (1 + 2*s2)) / sqrt(1 + 2*s2);

%% Report
fprintf('\n')
fprintf('  max |error|        : %.2e   (tol = %g)\n', max(abs(P-exact)), tol)
fprintf('  finest level       : L = %d  (%d cells)\n', info.L, nr*2^info.L)
fprintf('  MLMC cost          : %.3g\n', info.cost)
fprintf('  standard MC cost   : %.3g  (estimated, same tol, finest level)\n', ...
        info.mc_cost)
fprintf('  savings factor     : %.1f\n', info.mc_cost/info.cost)
fprintf('  wall-clock time    : %.2f s\n', elapsed)

figure('Name', 'Mean solution')
subplot(2,1,1)
plot(x, exp(-x.^2), 'k:', x, exact, 'b-', x, P, 'r.')
legend('u(x,0)', 'exact E[u(x,T)]', 'MLMC estimate', 'Location', 'NorthEast')
xlabel('x'); ylabel('u'); grid on
title('1D advection with random wave speed')
subplot(2,1,2)
plot(x, P - exact, 'r.-')
xlabel('x'); ylabel('error'); grid on

mlmc_plot(info, Nl, M, '1D advection with random wave speed')
