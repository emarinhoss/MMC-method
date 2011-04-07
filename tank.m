clear all; close all; clc;

tspan = [0 30];

gamma=4;
beta=2;

h0 = 1;

[t,h] = ode45('tankfill',tspan,h0,[],gamma,beta);
plot(t,h)