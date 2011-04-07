clear all; close all; clc;

s = ['g','m','b','r','k','y'];

%% Inputs
tfinal = 30;    % final time
dt = 0.1;      % timesteps to take
gamma = 4;
beta  = 2;
dts = 6;
%% Calculations
% interpolation times
tt = 0:dt:tfinal;
Plm1= zeros(1,length(tt));
Yl = Plm1;

for k=1:dts
    dtt = 2^(-k+1)*dt;
    n=round(tfinal/dtt);
    h = zeros(1,n+1);
    %%% Initial condition
    h(1) = 1;
    %%% -----------------
    t = 0;
    for i=2:length(h)
        t = t+dtt;
        h(i) = h(i-1) + dtt * (10+gamma*sin(t)-beta*sqrt(h(i-1)));
    end
    time = 0:dtt:tfinal;
    
    % calculate the Pl values for the current dt
    Pl = interp1(time,h,tt);
    
    % estimator summation
    Yl = Yl + (Pl - Plm1);
    
    % update P_{l-1}
    Plm1 = Pl;
    
    plot(time,h,s(k)), hold on
end

%Yl = Yl/dts;

plot(tt,Yl,'k')