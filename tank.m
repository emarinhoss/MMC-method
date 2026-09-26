%TANK  Deterministic solution of the tank-filling problem with ode45.
%
%   Solves  dh/dt = 10 + gamma*sin(t) - beta*sqrt(h),  h(0) = h0,
%   for fixed gamma and beta and plots h(t).  Useful as a reference for
%   the forward Euler solutions in DTANKFILL.
%
%   See also TANKFILL, DTANKFILL, MLMC_TANK.

clear; close all; clc;

tspan = [0 30];     % time interval
gamma = 4;          % amplitude of the periodic inflow
beta  = 2;          % outflow coefficient
h0    = 1;          % initial height

[t, h] = ode45(@(t,h) tankfill(t, h, gamma, beta), tspan, h0);

plot(t, h)
xlabel('t'); ylabel('h(t)'); grid on
title(sprintf('Tank filling, \\gamma = %g, \\beta = %g (ode45)', gamma, beta))
