function mlmc_plot(info, Nl, M, title_str)
%MLMC_PLOT  Standard MLMC diagnostic plots.
%
%   MLMC_PLOT(info, Nl, M, title_str)
%
%   Plots, against the level l,
%     1) log_M of the variance of P_l and of Y_l = P_l - P_{l-1},
%     2) log_M of |E[P_l]| and |E[Y_l]|,
%     3) the number of samples N_l used on each level.
%   For vector-valued outputs the maximum over the components is shown.
%
%   A healthy MLMC run shows V[Y_l] and |E[Y_l]| decaying linearly on the
%   log scale (the slopes are the strong and weak convergence rates) while
%   V[P_l] and |E[P_l]| stay roughly constant.
%
%   Inputs are the outputs of MLMC_ESTIMATE, the refinement factor M and a
%   title for the figure.
%
%   See also MLMC_ESTIMATE.

if nargin < 4, title_str = 'MLMC'; end

l    = 0:info.L;
logM = @(x) log(max(x, realmin)) / log(M);

figure('Name', title_str);

subplot(1,3,1)
plot(l, logM(max(info.VPl,[],2)), '*-', ...
     l(2:end), logM(max(info.Vl(2:end,:),[],2)), '*--')
xlabel('level l'); ylabel('log_M variance')
legend('P_l', 'P_l - P_{l-1}', 'Location', 'SouthWest')
grid on

subplot(1,3,2)
plot(l, logM(max(abs(info.EPl),[],2)), '*-', ...
     l(2:end), logM(max(abs(info.El(2:end,:)),[],2)), '*--')
xlabel('level l'); ylabel('log_M |mean|')
legend('P_l', 'P_l - P_{l-1}', 'Location', 'SouthWest')
grid on

subplot(1,3,3)
semilogy(l, Nl, '*-')
xlabel('level l'); ylabel('N_l')
grid on

if exist('sgtitle', 'file') || exist('sgtitle', 'builtin')
    sgtitle(title_str)
end
end
