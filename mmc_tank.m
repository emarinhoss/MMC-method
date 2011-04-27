clear all; close all; clc;
% initial number of temporal level L + 1
% L ==> initialy only 1 temporal value
L = -1;
% factor by which the timestep is refined
M = 4;
% Total integration interval
tfinal = 30;
% convergence tolerance
tol = 10e-4;
% 
N = 1e4;

% run until system converges
converged = 0;
while ~converged
    %
    % Step 1: define temporal levels for calculations
    %
    
    L = L+1; % increase L by 1 until system converges
    
    % timesteps
    hl = tfinal./(M.^(0:L));
    
    %
    % Step 2: Estimate VL using an initial
    %
    
    sums = MMC_tank(M,L,1e4);
    suml(1,L+1) = N;
    suml(2,L+1) = sums(1);
    suml(3,L+1) = sums(2);
    
    %
    % Step 3: Define optimal Nl
    %
    
    % Varinace Vl
    % V(X) = E(X^2) - (E(X))^2
    Vl = suml(3,:)./suml(1,:) - (suml(2,:)./suml(1,:)).^2;
    %  Optimal value of Nl
    Nl = ceil( 2*tol^(-2) * sqrt(Vl.*hl) * sum(sqrt(Vl.*hl)));
    
    %
    % Step 4: evaluate extra samples at each level
    %         as needed from new Nl
    
    Sums = MMC_tank(M,L,1e4);
    for l=0:L
        dNl = Nl(l+1)-suml(1,l+1);
        if dNl>0
            sums = MMC_tank(M,l,dNl);
            suml(1,l+1) = suml(1,l+1) + dNl;
            suml(2,l+1) = suml(2,l+1) + sums(1);
            suml(3,l+1) = suml(3,l+1) + sums(2);
        end
    end
    
    %
    % Step 5: Test for convergence
    %
    
    range = -1:0;
    if L>1 && M^L>=16
        con = M.^range.*suml(2,L+1+range)./suml(1,L+1+range);
        converged = (max(abs(con)) < (M-1)*tol/sqrt(2)) || (M^L>1024) ;
    end
    
end % while not converged loop

%
% evaluate multi-timestep estimator
%
P = sum( suml(2,1:L+1)./suml(1,1:L+1) );