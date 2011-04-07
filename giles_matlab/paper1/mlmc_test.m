%
% this code tests the multi-timestep Monte Carlo method
% on geometric Brownian motion with various payoffs
% and the Heston model
%

function mlmc_test

close all; clear all;

global option

for option = 1:5

  if option==1
    stitle = 'GBM model; European call';
  elseif option==2
    stitle = 'GBM model: Asian option';
  elseif option==3
    stitle = 'GBM model: lookback option';
  elseif option==4
    stitle = 'GBM model: digital call';
  elseif option==5
    stitle = 'Heston model: European call';
  end
  disp(stitle)

  if option==4
    N = 4000000;
  else
    N = 2000000;
  end

  M = 4;
  randn('state',0);

%
% first, convergence tests
%

  if option==3
    L = 0:5;
  else
    L = 0:4;
  end

  del1 = [];
  del2 = [];
  var1 = [];
  var2 = [];

  for l = L
    disp(sprintf('l = %d',l))
    sums = mlmc_l(M,l,N);
    del1 = [del1 sums(3)/N ];
    del2 = [del2 sums(1)/N ];
    var1 = [var1 sums(4)/N-(sums(3)/N)^2 ];
    var2 = [var2 sums(2)/N-(sums(1)/N)^2 ];
  end

  disp(sprintf('value = %f',del1(length(L))))

%
% second, mlmc complexity tests
%

  if option==3
    Eps = [ 0.002 0.001 0.0005 0.0002 0.0001 ];
  elseif option==4
    Eps = [ 0.005 0.002 0.001 0.0005 0.0002 ];
  else
    Eps = [ 0.001 0.0005 0.0002 0.0001 0.00005 ];
  end

  maxl = 0;

  for extrap = 0:1
    for i = 1:length(Eps)
      eps = Eps(i);
      [P, Nl] = mlmc(M,eps,@mlmc_l,extrap);
      l = length(Nl)-1;
      maxl = max(maxl,l);
      mlmc_cost(i,extrap+1) = (1+1/M)*sum(Nl.*M.^(0:l));
      std_cost(i,extrap+1)  = sum((2*var1((0:l)+1)/eps^2).*M.^(0:l));

      disp(sprintf('mlmc_cost = %d, std_cost = %d, savings = %f\n',...
         mlmc_cost(i,extrap+1),std_cost(i,extrap+1), ...
         std_cost(i,extrap+1)/ mlmc_cost(i,extrap+1)))
      Nls{i,extrap+1} = Nl;
    end
  end

  maxl = max(maxl,4);

%
% third, MSE evaluation
%

  if 1 & (option==1 | option==4 | option==5)
    if option==1
      P_exact = european_call(0.05,0.2,1,1,'value');
    elseif option==4
      P_exact = digital_call(0.05,0.2,1,1,'value');
    elseif option==5
      P_exact = heston_call(0.05,5,0.04,0.25,-0.5,1,1,0.04,1);
    end

    for extrap = 0:1
      for i = 1:length(Eps)
        eps = Eps(i);
        for k = 1:10
          [Pk(k), Nl] = mlmc(M,eps,@mlmc_l,extrap);
        end
        RMSE = sqrt( sum((Pk-P_exact).^2) / length(Pk) );
        disp(sprintf('RMSE, eps, ratio = %f %f %f',RMSE,eps,RMSE/eps))
      end
      disp(' ')
    end

  end

%
% plot figures
%

  nvert = 2;
  figure; pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.75*nvert]; set(gcf,'pos',pos);

  set(0,'DefaultAxesColorOrder',[0 0 0]);
  set(0,'DefaultAxesLineStyleOrder','-*|--*|:*|-.*')

  subplot(nvert,2,1)
  plot(L,log(var1)/log(M),L(2:end),log(var2(2:end))/log(M))
  xlabel('l'); ylabel('log_M variance'); %title(stitle)
  legend('P_l','P_l- P_{l-1}',3)
  axis([0 maxl -10 0])

  subplot(nvert,2,2)
  plot(L,log(abs(del1))/log(M), L(2:end),log(abs(del2(2:end)))/log(M), ...
                L(3:end),log(abs(del2(3:end)-del2(2:end-1)/M))/log(M))
  xlabel('l'); ylabel('log_M |mean|'); %title(stitle)
  legend('P_l','P_l- P_{l-1}','Y_l-Y_{l-1}/M',3)
  axis([0 maxl -12 0])

  set(0,'DefaultAxesLineStyleOrder','--o|--x|--d|--*|--s|:o|:x|:d|:*|:s');

  if nvert==1
    if option==1
      print('-deps2c','mlmc_gbm1a.eps')
    elseif option==2
      print('-deps2c','mlmc_gbm2a.eps')
    elseif option==3
      print('-deps2c','mlmc_gbm3a.eps')
    elseif option==4
      print('-deps2c','mlmc_gbm4a.eps')
    elseif option==5
      print('-deps2c','mlmc_hestona.eps')
    end
    figure; pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.75]; set(gcf,'pos',pos);
  end

  subplot(nvert,2,2*nvert-1)
  semilogy(0:length(Nls{5,1})-1,Nls{5,1}, ...
           0:length(Nls{4,1})-1,Nls{4,1}, ...
           0:length(Nls{3,1})-1,Nls{3,1}, ...
           0:length(Nls{2,1})-1,Nls{2,1}, ...
           0:length(Nls{1,1})-1,Nls{1,1}, ...
           0:length(Nls{5,2})-1,Nls{5,2}, ...
           0:length(Nls{4,2})-1,Nls{4,2}, ...
           0:length(Nls{3,2})-1,Nls{3,2}, ...
           0:length(Nls{2,2})-1,Nls{2,2}, ...
           0:length(Nls{1,2})-1,Nls{1,2});
  xlabel('l'); ylabel('N_l'); %title(stitle)
  if option==3
    legend('\epsilon=0.0001','\epsilon=0.0002','\epsilon=0.0005','\epsilon=0.001','\epsilon=0.002',1)
  elseif option==4
    legend('\epsilon=0.0002','\epsilon=0.0005','\epsilon=0.001','\epsilon=0.002','\epsilon=0.005',1)
  else
    legend('\epsilon=0.00005','\epsilon=0.0001','\epsilon=0.0002','\epsilon=0.0005','\epsilon=0.001',1)
  end

  if option==4
    axis([0 maxl 2e3 1e11])
  else
    axis([0 maxl 5e2 1e10])
  end

  set(0,'DefaultAxesLineStyleOrder','-*|-.*|--*|:*')

  subplot(nvert,2,2*nvert)
  loglog(Eps,Eps.^2.*std_cost(:,1)', ...
         Eps,Eps.^2.*std_cost(:,2)', ...
         Eps,Eps.^2.*mlmc_cost(:,1)', ...
         Eps,Eps.^2.*mlmc_cost(:,2)')
  xlabel('\epsilon'); ylabel('\epsilon^2 Cost'); %title(stitle)
    legend('Std MC','Std MC ext','MLMC','MLMC ext',3)
  if option==3
    axis([1e-4 2e-3 0.005 100])
  elseif option==4
    axis([2e-4 5e-3 0.2 200])
  else
    axis([5e-5 1e-3 0.005 10])
  end

  if nvert==1
    if option==1
      print('-deps2c','mlmc_gbm1b.eps')
    elseif option==2
      print('-deps2c','mlmc_gbm2b.eps')
    elseif option==3
      print('-deps2c','mlmc_gbm3b.eps')
    elseif option==4
      print('-deps2c','mlmc_gbm4b.eps')
    elseif option==5
      print('-deps2c','mlmc_hestonb.eps')
    end

  else
    if option==1
      print('-deps2c','mlmc_gbm1.eps')
    elseif option==2
      print('-deps2c','mlmc_gbm2.eps')
    elseif option==3
      print('-deps2c','mlmc_gbm3.eps')
    elseif option==4
      print('-deps2c','mlmc_gbm4.eps')
    elseif option==5
      print('-deps2c','mlmc_heston.eps')
    end
  end

end


%-------------------------------------------------------
%
% level l estimator
%

function sums = mlmc_l(M,l,N)

global option

T   = 1;
r   = 0.05;
sig = 0.2;

nf = M^l;
nc = nf/M;

hf = T/nf;
hc = T/nc;

sums(1:4) = 0;

for N1 = 1:10000:N
  N2 = min(10000,N-N1+1);

%
% GBM model
%
  if option<5
    X0 = 1;

    Xf = X0*ones(1,N2);
    Xc = Xf;

    Af = 0.5*hf*Xf;
    Ac = 0.5*hc*Xc;

    Mf = Xf;
    Mc = Xc;

    if l==0
      dWf = sqrt(hf)*randn(1,N2);
      Xf  = Xf + r*Xf*hf + sig*Xf.*dWf;
      Af = Af + 0.5*hf*Xf;
      Mf = min(Mf,Xf);
    else
      for n = 1:nc
        dWc = zeros(1,N2);
        for m = 1:M
          dWf = sqrt(hf)*randn(1,N2);
          dWc = dWc + dWf;
          Xf  = Xf + r*Xf*hf + sig*Xf.*dWf;
          Af  = Af + hf*Xf;
          Mf  = min(Mf,Xf);
        end
        Xc = Xc + r*Xc*hc + sig*Xc.*dWc;
        Ac = Ac + hc*Xc;
        Mc = min(Mc,Xc);
      end
      Af = Af - 0.5*hf*Xf;
      Ac = Ac - 0.5*hc*Xc;
    end

    if option==1
      Pf = max(0,Xf-1);
      Pc = max(0,Xc-1);
    elseif option==2
      Pf = max(0,Af-1);
      Pc = max(0,Ac-1);
    elseif option==3
      beta = 0.5826;  % special factor for offset correction
      Pf = Xf - Mf*(1-beta*sig*sqrt(hf));
      Pc = Xc - Mc*(1-beta*sig*sqrt(hc));
    elseif option==4
      Pf = 0.5*(sign(Xf-1)+1);
      Pc = 0.5*(sign(Xc-1)+1);
    end
%
% Heston model
%
  else
    X0 = [1; 0.04];
    Xf = X0*ones(1,N2);
    Xc = Xf;

    if l==0
      dWf = sqrt(hf)*randn(2,N2);
      Xf  = Xf + mu(Xf,hf)*hf + sig_dW(Xf,dWf,hf);

    else
      for n = 1:nc
        dWc = zeros(2,N2);
        for m = 1:M
          dWf = sqrt(hf)*randn(2,N2);
          dWc = dWc + dWf;
          Xf  = Xf + mu(Xf,hf)*hf + sig_dW(Xf,dWf,hf);
        end
        Xc = Xc + mu(Xc,hc)*hc + sig_dW(Xc,dWc,hc);
      end
    end

    Pf = max(0,Xf(1,:)-1);
    Pc = max(0,Xc(1,:)-1);
  end

  Pf = exp(-r*T)*Pf;
  Pc = exp(-r*T)*Pc;

  if l==0
    Pc=0;
  end

  sums(1) = sums(1) + sum(Pf-Pc);
  sums(2) = sums(2) + sum((Pf-Pc).^2);
  sums(3) = sums(3) + sum(Pf);
  sums(4) = sums(4) + sum(Pf.^2);
end

%--------------------

function m=mu(x,h)

%m = [ 0.05*x(1,:); ...
%       5*(0.04-x(2,:)) ];

m = [ 0.05*x(1,:); ...
       ((1-exp(-5*h))/h)*(0.04-x(2,:)) ];

%--------------------

function sigdW=sig_dW(x,dW,h)

dW(2,:) = -0.5*dW(1,:) + sqrt(0.75)*dW(2,:);

%sigdW = [ sqrt(max(0,x(2,:))).*x(1,:).*dW(1,:);  ...
%          0.25*sqrt(max(0,x(2,:))).*dW(2,:) ];

sigdW = [ sqrt(max(0,x(2,:))).*x(1,:).*dW(1,:);  ...
          exp(-5*h)*0.25*sqrt(max(0,x(2,:))).*dW(2,:) ];

