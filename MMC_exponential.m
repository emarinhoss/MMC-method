clear all; close all; clc;

y0 = 10;
tf = 1;
k = 3;
ks= 0.5;

% convergence tolerance
tol = 1e-2;

%% exact solution
exact = y0*exp(k*tf);

%% Monte Carlo
N =1e5;
dt=1e-3;
N2=0; sum_y=0; ssum_y=0;
% while 1
%     y = y0*ones(N,1);
%     r = k + ks.*randn(N,1);
%     t=0;
%     while t <= tf
%         y = y + dt*(r.*y);
%         t = t + dt;
%     end
% 
%     %ymean = mean(y);
%     clc
%     disp('Monte Carlo')
%     sum_y = sum_y + sum(y);
%     ssum_y= ssum_y + sum(y.*y);
%     N2 = N2 + N
%     ymean = sum_y / N2;
%     test_conv = abs(1-ymean/exact)
%     if (test_conv <= tol)
%         break
%     end
% end
% 
% yvari = ssum_y/N2 - ymean.*ymean;

%% Multilevel Monte Carlo
disp('Multilevel Monte Carlo')
% initial number of temporal level L + 1
% L ==> initialy only 1 temporal value
L = -1;
% factor by which the timestep is refined
M = 3; 

% run until system converges
converged = 0;
while ~converged
    %
    % Step 1: define temporal levels for calculations

    L = L+1 % increase L by 1 until system converges
    
    % timesteps
    hl = tf./(M.^(0:L));
    
    %
    % Step 2: Estimate VL using an initial
    
    sums = MMC_exp(M,L,N,tf,y0,k,ks);
    suml(1,L+1) = N;
    suml(2,L+1) = sums(1);
    suml(3,L+1) = sums(2);
    suml(4,L+1) = sums(3);
    suml(5,L+1) = sums(4);
    
    %
    % Step 3: Define optimal Nl
    
    % Varinace Vl
    % V(X) = E(X^2) - (E(X))^2
    Vl = abs(suml(3,:)./suml(1,:) - (suml(2,:)./suml(1,:)).^2);
    %  Optimal value of Nl
    Nl = ceil( 2*tol^(-2) * sqrt(Vl.*hl) * sum(sqrt(Vl.*hl)));
    
    %
    % Step 4: evaluate extra samples at each level
    %         as needed from new Nl
    
    for l=0:L
        dNl = Nl(l+1)-suml(1,l+1);
        if dNl>0
            sums = MMC_exp(M,l,dNl,tf,y0,k,ks);
            suml(1,l+1) = suml(1,l+1) + dNl;
            suml(2,l+1) = suml(2,l+1) + sums(1);
            suml(3,l+1) = suml(3,l+1) + sums(2);
            suml(4,l+1) = suml(4,l+1) + sums(3);
            suml(5,l+1) = suml(5,l+1) + sums(4);
        end
    end
    
    %
    % Step 5: Test for convergence
    % suml = abs(suml);
    range = -1:0;
    if L>1 && M^L>=16
        YL = M.^range.*suml(2,L+1+range)./suml(1,L+1+range);
        converged = (max(abs(YL)) < (M-1)*tol/sqrt(2)) || (M^L>1024) ;
    end
    
end % while not converged loop
%
% evaluate multi-timestep estimator
%
P = sum(suml(2,1:L+1)./suml(1,1:L+1));

% Multilevel Cost
C = Nl(1);
for k=2:length(Nl)
    l = k-1;
    C = C + Nl(k)*(M^l+M^(l-1));
end