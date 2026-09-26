function run_tests
%RUN_TESTS  Quick self-checks for the MLMC code (MATLAB or GNU Octave).
%
%   Runs each level sampler on a few samples to check the interface and
%   then runs small MLMC estimates for the problems with a known exact
%   mean, checking the error against a loose multiple of the tolerance.
%   Takes a few seconds.  Usage:  run_tests

rng(1);
nfail = 0;

% --- interface of the level samplers: 4-by-Q sums, Y_0 = P_0 -------------
s = mlmc_exponential_level(3, 0, 10, 1, 10, 3, 0.5);
nfail = nfail + check('exponential level: size', isequal(size(s), [4 1]));
nfail = nfail + check('exponential level: Y_0 = P_0', abs(s(1)-s(3)) < 1e-12);

s = mlmc_tank_level(4, 2, 10, 1, 1, 4, 0, 1);
nfail = nfail + check('tank level: size', isequal(size(s), [4 1]));
nfail = nfail + check('tank level: real', isreal(s));

s = mlmc_wave_level(1, 10, 0.5, 5, 16, -1, 0.1, 0.9);
nfail = nfail + check('wave level: size', isequal(size(s), [4 17]));

% --- batching: sums over N samples must not depend on the batch size -----
rng(2); a = mlmc_tank_level(4, 1, 25000, 1, 1, 4, 0, 1);
nfail = nfail + check('tank level: all batches counted', ...
                      abs(a(3)/25000 - 13) < 1);

% --- MLMC against exact solutions ----------------------------------------
quiet = struct('verbose', false);

y0 = 10; T = 1; k = 3; ks = 0.5; M = 3; tol = 1;
f = @(l,N) mlmc_exponential_level(M, l, N, T, y0, k, ks);
opts = quiet; opts.N0 = 1e3; opts.Lmin = 3;
P = mlmc_estimate(f, M, tol, opts);
exact = y0*exp(k*T + 0.5*(ks*T)^2);
nfail = nfail + check(sprintf('exponential MLMC (err %.2g, tol %g)', ...
                      abs(P-exact), tol), abs(P-exact) < 3*tol);

T = 0.5; R = 5; nr = 32; cm = -1; cs = 0.2; tol = 5e-3;
f = @(l,N) mlmc_wave_level(l, N, T, R, nr, cm, cs, 0.9);
opts = quiet; opts.N0 = 100; opts.cost_exp = 2;
P = mlmc_estimate(f, 2, tol, opts);
x = linspace(-R, R, nr+1); s2 = (cs*T)^2;
exact = exp(-(x - cm*T).^2/(1 + 2*s2))/sqrt(1 + 2*s2);
err = max(abs(P - exact));
nfail = nfail + check(sprintf('wave MLMC (max err %.2g, tol %g)', err, tol), ...
                      err < 3*tol);

% --- convergence-test utility ---------------------------------------------
f = @(l,N) mlmc_tank_level(4, l, N, 1, 1, 4, 0, 1);
Vl = estimate_Vl(f, 1:3, 1e4);
nfail = nfail + check('estimate_Vl: variance decays', all(diff(Vl) < 0));

if nfail == 0
    fprintf('\nAll tests passed.\n')
else
    error('%d test(s) failed.', nfail)
end
end

function failed = check(name, ok)
%CHECK  Print a PASS/FAIL line; return 1 on failure.
if ok
    fprintf('PASS  %s\n', name)
else
    fprintf('FAIL  %s\n', name)
end
failed = ~ok;
end
