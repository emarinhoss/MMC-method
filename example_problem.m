%EXAMPLE_PROBLEM  Plain Monte Carlo integration of x^3 on [a, b].
%
%   Estimates  I = 1/(b-a) * integral_a^b x^3 dx = E[X^3],  X ~ U(a, b),
%   with standard Monte Carlo for an increasing number of samples N, and
%   compares with the exact value.  The error decays like
%   sqrt(Var[X^3]/N), the O(N^(-1/2)) rate that multilevel Monte Carlo
%   improves upon for problems that also require a discretisation.
%
%   For a = 0, b = 1:  E[X^3] = 1/4,  Var[X^3] = 1/7 - 1/16 = 9/112.

clear; close all; clc;
rng(0);                             % reproducible results

a    = 0;                           % integration interval [a, b]
b    = 1;
Nmax = 1e6;                         % largest sample size
NN   = round(logspace(2, log10(Nmax), 40));   % sample sizes to test

% exact mean and variance of X^3 for X ~ U(a,b)
Ex3   = (b^4 - a^4) / (4*(b-a));
Ex6   = (b^7 - a^7) / (7*(b-a));
exact = Ex3;
var_exact = Ex6 - Ex3^2;

mcarlo = zeros(size(NN));           % Monte Carlo estimates
mvar   = zeros(size(NN));           % sample variances
for i = 1:numel(NN)
    x = a + (b-a)*rand(NN(i), 1);
    f = x.^3;
    mcarlo(i) = mean(f);
    mvar(i)   = var(f);
end

subplot(2,1,1)
loglog(NN, abs(mcarlo - exact), 'o-', NN, sqrt(var_exact./NN), 'k--')
xlabel('N'); ylabel('|error|'); grid on
legend('Monte Carlo', 'sqrt(Var/N)', 'Location', 'SouthWest')
title('Monte Carlo estimate of E[X^3]')

subplot(2,1,2)
semilogx(NN, mvar, 'o-', NN, var_exact*ones(size(NN)), 'k--')
xlabel('N'); ylabel('variance'); grid on
legend('sample variance', 'exact variance', 'Location', 'NorthEast')
