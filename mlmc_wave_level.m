function sums = mlmc_wave_level(l, N, T, R, nr, c_mean, c_std, cfl)
%MLMC_WAVE_LEVEL  Level-l MLMC sampler for the 1D linear advection (one-way
%   wave) equation with a random wave speed.
%
%   sums = MLMC_WAVE_LEVEL(l, N, T, R, nr, c_mean, c_std, cfl)
%
%   Model problem
%       u_t + c*u_x = 0,   x in [-R, R],   u(x,0) = exp(-x^2),
%       c ~ Normal(c_mean, c_std^2),
%   with quantity of interest P = u(x,T) evaluated at the nr+1 points of
%   the coarsest grid.  The exact solution is u(x,t) = exp(-(x - c*t)^2).
%
%   Discretisation: first-order upwind finite differences with explicit
%   Euler time stepping.  Level l uses nr*2^l cells and n0*2^l time steps,
%   i.e. the mesh size and the time step are both halved from one level to
%   the next (refinement factor M = 2).  The upwind direction is chosen
%   per sample from the sign of c, and the solution at the inflow boundary
%   is held at its initial value (exp(-R^2), essentially zero).  The fine
%   and coarse solutions of each sample share the same c.
%
%   The time step on each level is fixed (independent of the sample), so
%   that all samples in a batch can be advanced together.  It is chosen so
%   that the CFL number |c|*dt/dx equals CFL for |c| = |c_mean| + 5*c_std;
%   faster samples (probability < 1e-6) have a slightly larger CFL number,
%   which is still stable as long as it stays below 1.
%
%   Inputs
%     l       level (0, 1, 2, ...)
%     N       number of samples
%     T       final time
%     R       half-width of the domain [-R, R]
%     nr      number of cells on the coarsest grid (level 0)
%     c_mean  mean of the random wave speed
%     c_std   standard deviation of the random wave speed
%     cfl     target CFL number (< 1)
%
%   Output
%     sums    4-by-(nr+1) array [sum(Y); sum(Y.^2); sum(P_l); sum(P_l.^2)]
%             at each output point, as required by MLMC_ESTIMATE.
%
%   See also MLMC_ESTIMATE, MLMC_WAVE.

c_max = abs(c_mean) + 5*c_std;
dx0   = 2*R/nr;
n0    = ceil(T*c_max/(cfl*dx0));    % number of time steps on level 0

batch = 1e4;     % samples are processed in batches to limit memory use
sums  = zeros(4, nr+1);

for N1 = 1:batch:N
    N2 = min(batch, N-N1+1);

    c  = c_mean + c_std*randn(N2,1);         % random wave speed per sample

    Pf = advect(c, T, R, nr*2^l, n0*2^l);    % fine solution
    Pf = Pf(:, 1:2^l:end);                   % restrict to the output points
    if l == 0
        Pc = zeros(size(Pf));                % P_{-1} = 0, so Y_0 = P_0
    else
        Pc = advect(c, T, R, nr*2^(l-1), n0*2^(l-1));
        Pc = Pc(:, 1:2^(l-1):end);
    end

    Y = Pf - Pc;
    sums(1,:) = sums(1,:) + sum(Y,    1);
    sums(2,:) = sums(2,:) + sum(Y.^2, 1);
    sums(3,:) = sums(3,:) + sum(Pf,   1);
    sums(4,:) = sums(4,:) + sum(Pf.^2,1);
end
end

% -------------------------------------------------------------------------
function u = advect(c, T, R, nx, nt)
%ADVECT  Solve u_t + c*u_x = 0 on [-R,R] up to time T with first-order
%   upwinding, nx cells and nt time steps.  c is a column vector with one
%   wave speed per sample; u is numel(c)-by-(nx+1) (one row per sample).
x   = linspace(-R, R, nx+1);
dx  = 2*R/nx;
dt  = T/nt;
lam = c*dt/dx;                        % signed CFL number of each sample
neg = lam < 0;                        % c < 0: information travels left
pos = ~neg;

u = repmat(exp(-x.^2), numel(c), 1);  % initial condition
for n = 1:nt
    du = diff(u, 1, 2);               % u(j+1) - u(j), size N-by-nx
    % c < 0: forward difference, the right boundary is the inflow
    u(neg, 1:end-1) = u(neg, 1:end-1) - lam(neg).*du(neg,:);
    % c > 0: backward difference, the left boundary is the inflow
    u(pos, 2:end)   = u(pos, 2:end)   - lam(pos).*du(pos,:);
    % inflow boundary nodes are left unchanged
end
end
