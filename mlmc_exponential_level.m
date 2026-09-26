function sums = mlmc_exponential_level(M, l, N, T, y0, k_mean, k_std)
%MLMC_EXPONENTIAL_LEVEL  Level-l MLMC sampler for exponential growth with
%   a random rate.
%
%   sums = MLMC_EXPONENTIAL_LEVEL(M, l, N, T, y0, k_mean, k_std)
%
%   Model problem
%       dy/dt = k*y,   y(0) = y0,   k ~ Normal(k_mean, k_std^2),
%   with quantity of interest P = y(T).  The ODE is integrated with the
%   forward Euler method using M^l steps on level l.
%
%   The fine (M^l steps) and coarse (M^(l-1) steps) solutions of each
%   sample share the same random rate k, which is what makes the variance
%   of the correction Y_l = P_l - P_{l-1} small.
%
%   Inputs
%     M       refinement factor between levels
%     l       level (0, 1, 2, ...)
%     N       number of samples
%     T       final time
%     y0      initial condition
%     k_mean  mean of the random growth rate
%     k_std   standard deviation of the random growth rate
%
%   Output
%     sums    4-by-1 vector [sum(Y); sum(Y.^2); sum(P_l); sum(P_l.^2)],
%             as required by MLMC_ESTIMATE.
%
%   See also MLMC_ESTIMATE, MLMC_EXPONENTIAL.

nf = M^l;        % number of fine time steps
nc = nf/M;       % number of coarse time steps (only used for l > 0)
hf = T/nf;       % fine time step
hc = T/nc;       % coarse time step

batch = 1e5;     % samples are processed in batches to limit memory use
sums = zeros(4,1);

for N1 = 1:batch:N
    N2 = min(batch, N-N1+1);

    k  = k_mean + k_std*randn(1,N2);    % random growth rate, one per sample
    Pf = y0*ones(1,N2);                 % fine solution
    Pc = Pf;                            % coarse solution

    if l == 0
        Pf = Pf + hf*(k.*Pf);           % a single Euler step
        Pc = zeros(1,N2);               % P_{-1} = 0, so Y_0 = P_0
    else
        for n = 1:nc
            % M fine steps for every coarse step
            for m = 1:M
                Pf = Pf + hf*(k.*Pf);
            end
            Pc = Pc + hc*(k.*Pc);
        end
    end

    Y = Pf - Pc;
    sums(1) = sums(1) + sum(Y);
    sums(2) = sums(2) + sum(Y.^2);
    sums(3) = sums(3) + sum(Pf);
    sums(4) = sums(4) + sum(Pf.^2);
end
end
