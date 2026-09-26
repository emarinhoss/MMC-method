# MMC-method: Multilevel Monte Carlo in MATLAB

MATLAB / GNU Octave implementations of the **Multilevel Monte Carlo (MLMC)**
method of M. B. Giles, applied to a few model problems with random
coefficients: an ODE with a random growth rate, a tank-filling ODE with a
random outflow, and a 1D advection (one-way wave) PDE with a random wave
speed.

## Background

We want $\mathbb{E}[P]$, where $P$ is a quantity of interest (e.g. the
solution at the final time) that can only be computed approximately, by a
discretisation with step size $h_l = h_0 M^{-l}$ on level $l$. Standard
Monte Carlo on the finest level $L$ needs $O(\varepsilon^{-2})$ samples,
each of which is expensive. MLMC rewrites the expectation as a telescoping
sum

$$
\mathbb{E}[P_L] = \mathbb{E}[P_0] + \sum_{l=1}^{L} \mathbb{E}[P_l - P_{l-1}],
$$

and estimates each term independently. Because the fine and coarse
solutions $P_l$ and $P_{l-1}$ of one sample use the **same random input**,
the correction $Y_l = P_l - P_{l-1}$ has a small variance $V_l$, so only a
few samples are needed on the expensive fine levels. Most of the work is
done on the cheap coarse levels.

The adaptive algorithm (`mlmc_estimate.m`, after Giles 2008) is:

1. Add a level $L$ and take $N_0$ initial samples on it.
2. Estimate the variance $V_l$ of $Y_l$ on every level.
3. Choose the number of samples per level that minimises the cost for a
   sampling variance of $\varepsilon^2/2$:
   $N_l = \left\lceil 2\varepsilon^{-2}\sqrt{V_l/C_l}\,\sum_k \sqrt{V_k C_k}\right\rceil$,
   where $C_l$ is the cost of one sample on level $l$.
4. Take the extra samples needed on each level.
5. Stop when the estimated bias, based on
   $\max\left(|\mathbb{E}[Y_{L-1}]|/M,\ |\mathbb{E}[Y_L]|\right)$, is below
   $(M-1)\varepsilon/\sqrt{2}$. Otherwise go back to step 1.

The result $\sum_l \overline{Y_l}$ then has a root-mean-square error of
about $\varepsilon$.

## Requirements

* MATLAB R2016b or newer (the code uses implicit expansion), **or**
* a recent GNU Octave (tested with Octave 8.4).

You don't need any toolboxes.

## Quick start

From the repository root, in MATLAB or Octave:

```matlab
run_tests            % a few seconds of self-checks
mlmc_exponential     % ODE with a random growth rate
mlmc_tank            % tank-filling ODE with a random outflow
mlmc_wave            % 1D advection PDE with a random wave speed
```

Each driver prints the MLMC estimate, the error against an exact or
reference value, and the cost compared with standard Monte Carlo. It then
draws the usual MLMC diagnostics: variance and mean of $P_l$ and $Y_l$
against the level, and $N_l$ against the level.

## Repository layout

| File | Description |
|------|-------------|
| `mlmc_estimate.m` | Generic adaptive MLMC driver: takes a level-sampler function handle and returns the estimate, $N_l$ and per-level diagnostics. |
| `estimate_Vl.m` | Convergence test: mean and variance of $Y_l$ and $P_l$ for a fixed number of samples on a range of levels. Use it to check the variance decay before running MLMC. |
| `mlmc_plot.m` | Standard MLMC diagnostic plots. |
| `mlmc_exponential.m` | Driver: $y' = k y$, $k \sim N(3, 0.5^2)$. Compared with the exact mean $y_0 e^{\mu T + \sigma^2 T^2/2}$. |
| `mlmc_exponential_level.m` | Level sampler for the problem above (forward Euler, $M^l$ steps). |
| `mlmc_tank.m` | Driver: tank filling $h' = 10 + \gamma\sin t - \beta\sqrt{h}$, $\beta \sim N(0,1)$. Compared with a Gauss–Hermite / `ode45` reference. |
| `mlmc_tank_level.m` | Level sampler for the tank problem (forward Euler, $M^l$ steps). |
| `mlmc_wave.m` | Driver: $u_t + c u_x = 0$, $u(x,0) = e^{-x^2}$, $c \sim N(-1, 0.2^2)$. Compared with the exact mean solution. |
| `mlmc_wave_level.m` | Level sampler for the advection problem (first-order upwind, mesh and time step halved per level). |
| `tank.m`, `tankfill.m` | Deterministic tank problem solved with `ode45` (`tankfill.m` is the right-hand side). |
| `dtankfill.m` | Deterministic illustration of the telescoping sum: forward Euler with successively halved time steps. |
| `example_problem.m` | Plain Monte Carlo integration of $x^3$ on $[0,1]$, showing the $N^{-1/2}$ error decay. |
| `run_tests.m` | Quick self-checks of all samplers and of MLMC against exact solutions. |
| `giles_matlab/` | Reference MATLAB codes by M. B. Giles accompanying his MLMC papers, kept unmodified (see below). |

## Adding your own problem

Write a level sampler with the signature

```matlab
sums = my_level(l, N)
```

It draws `N` independent samples on level `l` and returns a `4-by-Q` array:

| Row | Contents |
|-----|----------|
| `sums(1,:)` | $\sum Y_l$, where $Y_l = P_l - P_{l-1}$ and $Y_0 = P_0$ |
| `sums(2,:)` | $\sum Y_l^2$ |
| `sums(3,:)` | $\sum P_l$ |
| `sums(4,:)` | $\sum P_l^2$ |

`Q` is the number of outputs: 1 for a scalar, or e.g. the number of grid
points for a PDE solution. Compute $P_l$ and $P_{l-1}$ with the **same**
random inputs. Then call

```matlab
opts = struct('N0', 1e3, 'Lmin', 2, 'Lmax', 10, 'cost_exp', 1);
[P, Nl, info] = mlmc_estimate(@my_level, M, tol, opts);
mlmc_plot(info, Nl, M, 'My problem');
```

`cost_exp` (call it $c$) sets how the cost grows with the level,
$C_l = M^{c\,l}$. Use 1 when only the time step is refined, and 2 for a 1D PDE refined in
both space and time. See `help mlmc_estimate` for all options and outputs.

## Example results

Output of the three drivers with the default parameters and `rng(0)`
(Octave 8.4). The savings compare the MLMC cost $\sum_l N_l C_l$ with the
estimated cost $2\,\mathbb{V}[P_L]\,\varepsilon^{-2} C_L$ of standard Monte
Carlo on the finest level:

| Problem | tol (RMS) | Error | Finest level | Estimated MLMC savings vs. standard MC |
|---------|-----------|-------|--------------|------------------------------------|
| Exponential growth, `mlmc_exponential` | 0.5 | 0.26 | $L=8$ ($3^8$ steps) | ≈ 58× |
| Tank filling, `mlmc_tank` | 1e-3 | 1.1e-3 | $L=6$ ($4^6$ steps) | ≈ 790× |
| 1D advection, `mlmc_wave` | 2e-3 | 1.5e-3 (max over $x$) | $L=6$ (4096 cells) | ≈ 110× |

Because `tol` bounds the *root-mean-square* error, a single run can
occasionally land slightly above it, as in the tank run.

## Giles' reference codes

`giles_matlab/` contains the original MATLAB programs by Mike Giles:

* `paper1/`: *Multilevel Monte Carlo path simulation*: `mlmc.m` (driver),
  `mlmc_test.m` (GBM with European, Asian, lookback and digital options,
  and the Heston model), analytic option prices, and the paper's figures
  (`.eps`).
* `paper2/`: *Improved multilevel Monte Carlo convergence using the
  Milstein scheme*: `mlmc.m` and `mlmc_test2.m` (the same GBM payoffs plus a
  barrier option, and the Heston model, with the Milstein scheme), plus
  the paper's figures.

The code in the repository root follows the same algorithm and the same
`sums` convention. It uses function handles instead of globals, and it
supports vector-valued outputs.

## References

1. M. B. Giles, "Multilevel Monte Carlo path simulation",
   *Operations Research* 56(3):607–617, 2008.
   [doi:10.1287/opre.1070.0496](https://doi.org/10.1287/opre.1070.0496)
2. M. B. Giles, "Improved multilevel Monte Carlo convergence using the
   Milstein scheme", in *Monte Carlo and Quasi-Monte Carlo Methods 2006*,
   pp. 343–358, Springer.
3. M. B. Giles, "Multilevel Monte Carlo methods", *Acta Numerica*
   24:259–328, 2015.
   [doi:10.1017/S096249291500001X](https://doi.org/10.1017/S096249291500001X)
