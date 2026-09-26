function dhdt = tankfill(t, h, gamma, beta)
%TANKFILL  Right-hand side of the tank-filling ODE.
%
%   dhdt = TANKFILL(t, h, gamma, beta) returns
%       dh/dt = 10 + gamma*sin(t) - beta*sqrt(h)
%   for the height h of liquid in a tank with a periodic inflow of
%   amplitude gamma and an outflow proportional to sqrt(h) (Torricelli's
%   law) with coefficient beta.  Suitable for ode45 via an anonymous
%   function, e.g.
%
%       [t, h] = ode45(@(t,h) tankfill(t, h, 4, 2), [0 30], 1);
%
%   See also TANK, MLMC_TANK_LEVEL.

dhdt = 10 + gamma*sin(t) - beta*sqrt(h);
end
