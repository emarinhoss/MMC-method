%DTANKFILL  Telescoping-sum illustration on the deterministic tank problem.
%
%   Solves  dh/dt = 10 + gamma*sin(t) - beta*sqrt(h),  h(0) = 1,
%   with forward Euler for a sequence of time steps dt, dt/2, dt/4, ...
%   and accumulates the corrections P_l - P_{l-1} on a common time grid.
%   By construction the telescoping sum
%       P_0 + (P_1 - P_0) + ... + (P_L - P_{L-1}) = P_L
%   reproduces the finest solution; the shrinking size of the corrections
%   is what multilevel Monte Carlo exploits (few samples are needed for the
%   small, expensive corrections).
%
%   See also TANK, TANKFILL, MLMC_TANK.

clear; close all; clc;

%% Inputs
tfinal = 30;        % final time
dt     = 0.1;       % coarsest time step (also the output grid spacing)
gamma  = 4;         % amplitude of the periodic inflow
beta   = 2;         % outflow coefficient
nlev   = 6;         % number of levels (time steps dt, dt/2, ..., dt/2^5)

colors = lines(nlev);

%% Calculations
tt   = 0:dt:tfinal;             % common output grid
Plm1 = zeros(size(tt));         % P_{l-1} on the output grid (P_{-1} = 0)
Yl   = zeros(size(tt));         % running telescoping sum
corr = zeros(1, nlev);          % max |P_l - P_{l-1}| for each level

figure; hold on
for k = 1:nlev
    dtk  = dt/2^(k-1);          % time step on this level
    n    = round(tfinal/dtk);   % number of steps
    time = (0:n)*dtk;

    % forward Euler
    h    = zeros(1, n+1);
    h(1) = 1;                   % initial condition
    for i = 1:n
        h(i+1) = h(i) + dtk*(10 + gamma*sin(time(i)) - beta*sqrt(h(i)));
    end

    % P_l on the common output grid (every 2^(k-1)-th point)
    Pl = h(1:2^(k-1):end);

    % accumulate the correction P_l - P_{l-1}
    corr(k) = max(abs(Pl - Plm1));
    Yl   = Yl + (Pl - Plm1);
    Plm1 = Pl;

    plot(time, h, 'Color', colors(k,:), ...
         'DisplayName', sprintf('dt = %g', dtk))
end

plot(tt, Yl, 'k--', 'DisplayName', 'telescoping sum')
xlabel('t'); ylabel('h(t)'); grid on; legend('Location', 'SouthEast')
title('Forward Euler with successively halved time steps')

fprintf('level   dt          max|P_l - P_{l-1}|   (P_{-1} = 0)\n')
for k = 1:nlev
    fprintf('%3d     %-10g  %.3e\n', k-1, dt/2^(k-1), corr(k))
end
