function sums = mlmc_tank_level(M, l, N, T, h0, gamma, beta_mean, beta_std)
%MLMC_TANK_LEVEL  Level-l MLMC sampler for the tank-filling problem with a
%   random outflow coefficient.
%
%   sums = MLMC_TANK_LEVEL(M, l, N, T, h0, gamma, beta_mean, beta_std)
%
%   Model problem (height h of liquid in a tank)
%       dh/dt = 10 + gamma*sin(t) - beta*sqrt(h),   h(0) = h0,
%       beta ~ Normal(beta_mean, beta_std^2),
%   with quantity of interest P = h(T).  The ODE is integrated with the
%   forward Euler method using M^l steps on level l.  The fine and coarse
%   paths of each sample share the same beta.
%
%   Inputs
%     M          refinement factor between levels
%     l          level (0, 1, 2, ...)
%     N          number of samples
%     T          final time
%     h0         initial height
%     gamma      amplitude of the periodic inflow
%     beta_mean  mean of the random outflow coefficient
%     beta_std   standard deviation of the random outflow coefficient
%
%   Output
%     sums       4-by-1 vector [sum(Y); sum(Y.^2); sum(P_l); sum(P_l.^2)],
%                as required by MLMC_ESTIMATE.
%
%   Note: forward Euler is only stable here while dt*beta/(2*sqrt(h)) < 2,
%   so very coarse levels on long time intervals can blow up.  The height
%   is clipped at zero inside sqrt() so that an overshoot cannot produce
%   complex numbers.
%
%   See also MLMC_ESTIMATE, MLMC_TANK, TANKFILL.

nf = M^l;        % number of fine time steps
nc = nf/M;       % number of coarse time steps (only used for l > 0)
hf = T/nf;       % fine time step
hc = T/nc;       % coarse time step

% right-hand side of the ODE, vectorised over samples
rhs = @(t, h, beta) 10 + gamma*sin(t) - beta.*sqrt(max(h,0));

batch = 1e4;     % samples are processed in batches to limit memory use
sums = zeros(4,1);

for N1 = 1:batch:N
    N2 = min(batch, N-N1+1);

    beta = beta_mean + beta_std*randn(1,N2);
    Pf = h0*ones(1,N2);          % fine solution
    Pc = Pf;                     % coarse solution
    tf = 0;                      % fine time
    tc = 0;                      % coarse time

    if l == 0
        Pf = Pf + hf*rhs(tf, Pf, beta);
        Pc = zeros(1,N2);        % P_{-1} = 0, so Y_0 = P_0
    else
        for n = 1:nc
            % M fine steps for every coarse step
            for m = 1:M
                Pf = Pf + hf*rhs(tf, Pf, beta);
                tf = tf + hf;
            end
            Pc = Pc + hc*rhs(tc, Pc, beta);
            tc = tc + hc;
        end
    end

    Y = Pf - Pc;
    sums(1) = sums(1) + sum(Y);
    sums(2) = sums(2) + sum(Y.^2);
    sums(3) = sums(3) + sum(Pf);
    sums(4) = sums(4) + sum(Pf.^2);
end
end
