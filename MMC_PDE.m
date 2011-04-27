%% Multilevel Monte Carlo
clear all; close all; clc;

N =1e4;     % initial sample size
tol = 1e-1;    % tolerance
R = 5;      % half domain
n=64;       % coarses resolution
T=0.1;      % final time

disp('Multilevel Monte Carlo')
% initial number of temporal level L + 1
% L ==> initialy only 1 temporal value
L = -1;
% factor by which the timestep is refined
M = 2; 

suml = zeros(4,n+1,20);
nums = zeros(1,20);
% run until system converges
converged = 0;
while ~converged
    %
    % Step 1: define temporal levels for calculations

    L = L+1 % increase L by 1 until system converges
    
    % timesteps
    hl = T./(M.^(0:L));
    
    %
    % Step 2: Estimate VL using an initial
    
    sums = mmc_wave(M,L,N,T,R,n);
    nums(L+1) = N;
    suml(1,:,L+1) = sums(1,:);
    suml(2,:,L+1) = sums(2,:);
    suml(3,:,L+1) = sums(3,:);
    suml(4,:,L+1) = sums(4,:);
    
    %
    % Step 3: Define optimal Nl
    
    % Varinace Vl
    % V(X) = E(X^2) - (E(X))^2
    Vl = zeros(L+1,n+1);
    for m=1:L+1
        Vl(m,:) = Vl(m,:) + abs(suml(2,:,m)/nums(m) - (suml(1,:,m)/nums(m)).^2);
    end
    V = max(Vl');
    %  Optimal value of Nl
    Nl = ceil( 2*tol^(-2) * sqrt(V.*hl) * sum(sqrt(V.*hl)));
    
    %
    % Step 4: evaluate extra samples at each level
    %         as needed from new Nl
    
    for l=0:L
        dNl = Nl(l+1)-suml(1,l+1);
        if dNl>0
            sums = mmc_wave(M,l,dNl,T,R,n);
            nums(l+1)     = nums(l+1)     + dNl;
            suml(1,:,l+1) = suml(1,:,l+1) + sums(1,:);
            suml(2,:,l+1) = suml(2,:,l+1) + sums(2,:);
            suml(3,:,l+1) = suml(3,:,l+1) + sums(3,:);
            suml(4,:,l+1) = suml(4,:,l+1) + sums(4,:);
        end
    end
    
    %
    % Step 5: Test for convergence

    %range = -1:0;
    if L >= 1
        for k=-1:0
            YL(k+2,:) = M.^k.*suml(1,:,L+1+k)./nums(L+1+k);
        end
        converged = (max(max(abs(YL'))) < (M-1)*tol/sqrt(2) || (L>20));
    end
    
end % while not converged loop
%
% evaluate multi-timestep estimator
%
P = suml(1,:,1)/nums(1);
for m=2:L+1
    P = P + suml(1,:,m)/nums(m);
end