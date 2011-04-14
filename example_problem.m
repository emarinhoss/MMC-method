% Testcase to compare Monte Carlo method with
% Multi level and the exact solution.
%
% funtion: y = x^3
clear all; clc;
%

N = 1e6;
a = 0;
b = 1;

% used to store the values
l=1;
% sampling increment
inc = 1e3; NN = inc:inc:N;
exact = zeros(1,length(NN)); mcarlo = exact; evar=exact; mvar=exact;
for m=NN
    %% exact solution
    % exact integration interval
%     dx = (b-a)/m;
%     
%     v = linspace(a,b-0.5*dx,m);
%     for k=1:m
%         x = v(k)+0.5*dx;
%         exact(l) = exact(l) + x*x*x*dx;
%         evar(l) = evar(l) + x*x*x*x*x*x;
%     end
%     exact(l) = (b-a).^(-1)*exact(l);
%     evar(l) = evar(l)/m-exact(l)*exact(l);
    
    %% Monte-Carlo
    r = a + (b-a).*rand(m,1);
    mcarlo(l) = 1/m*sum(r.*r.*r);
    %mvar(l) = sum(r.*r.*r.*r.*r.*r)/m - mcarlo(l)*mcarlo(l);
    l = l+1;
end

subplot(2,1,1),plot(NN,(mcarlo-0.25))
title('Error')
% subplot(2,1,2),plot(NN,evar,NN,mvar)
% title('Variance')
% legend('Exact Solution','Monte-Carlo')