%
% This code tests the multilevel Monte Carlo method using
% the Milstein discretisation.

% The Geometric Brownian Motion results, with a variety
% of payoffs, appeared in 'Improved multilevel Monte Carlo
% convergence using the Milstein scheme'. pp.343-358 in 
% Monte Carlo and Quasi-Monte Carlo Methods 2006, Springer, 2007. 
% http://people.maths.ox.ac.uk/~gilesm/psfiles/mcqmc06.pdf
%
% The results for the Heston model have not been published.
%
% copyright Mike Giles, 2006/7
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
    stitle = 'GBM model: barrier option';
  elseif option==6
    stitle = 'Heston model: European call';
  end
  disp(stitle)

  if option==4 | option==5 | option==6
    N = 100000;
  else
    N = 10000;
  end

  M = 2;

  randn('state',0);
  rand('state',0);

%
% first, convergence tests
%

  if option==3
    L = 0:10;
  else
    L = 0:8;
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
    var1 = [var1 max(1e-10, sums(4)/N-(sums(3)/N)^2) ];
    var2 = [var2 max(1e-10, sums(2)/N-(sums(1)/N)^2) ];

%
% this is a check on whether equality (5) in the MCQMC06 paper
% is satisfied to within statistical error; the output is 
% difference in estimated expectations / standard deviation
%
    if (option==3 | option==4 | option==5) & l>0
      disp(sprintf('err/s.d. = %f',abs(sums(5)/N) / sqrt((sums(6)/N-(sums(5)/N)^2)/N)))
    end
  end

%
% second, mlmc complexity tests
%

  randn('state',0);
  rand('state',0);
  
  if option==4
    Eps = [ 0.002 0.001 0.0005 0.0002 0.0001 ];
  else
    Eps = [ 0.001 0.0005 0.0002 0.0001 0.00005];
  end

  for i = 1:length(Eps)
    eps = Eps(i);
    [P, Nl] = mlmc(M,eps,@mlmc_l);
    l = length(Nl)-1;
    mlmc_cost(i) = sum(Nl.*M.^(0:l));
    std_cost(i)  = ceil(2*var1(l+1)/eps^2) * M^l;
    disp(sprintf('mlmc_cost = %d, std_cost = %d, savings = %f\n',...
         mlmc_cost(i),std_cost(i),std_cost(i)/ mlmc_cost(i)))
    Nls{i} = Nl;
  end

  l = max(l,4);
  l = max(L);

%
% plot figures
%

  nvert = 2;
  figure; pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.75*nvert]; set(gcf,'pos',pos);

  set(0,'DefaultAxesColorOrder',[0 0 0]);
  set(0,'DefaultAxesLineStyleOrder','-*|--*|-.*|:*')

  subplot(nvert,2,1)
  plot(L,log(var1)/log(M),L(2:end),log(var2(2:end))/log(M))
  xlabel('l'); ylabel('log_2 variance'); %title(stitle)
  legend('P_l','P_l- P_{l-1}',3)
  axis([0 l -60/M -0])


  subplot(nvert,2,2)
  plot(L,log(abs(del1))/log(M), L(2:end),log(abs(del2(2:end)))/log(M))
  xlabel('l'); ylabel('log_2 |mean|'); %title(stitle)
  legend('P_l','P_l- P_{l-1}',3)
  axis([0 l -40/M 0])

  set(0,'DefaultAxesLineStyleOrder','-o|-x|-d|-*|-s');

  if nvert==1
    if option==1
      print('-deps2c','mlmc2_gbm1a.eps')
    elseif option==2
      print('-deps2c','mlmc2_gbm2a.eps')
    elseif option==3
      print('-deps2c','mlmc2_gbm3a.eps')
    elseif option==4
      print('-deps2c','mlmc2_gbm4a.eps')
    elseif option==5
      print('-deps2c','mlmc2_gbm5a.eps')
    elseif option==6
      print('-deps2c','mlmc2_hestona.eps')
    end
    figure; pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.75]; set(gcf,'pos',pos);
  end

  subplot(nvert,2,2*nvert-1)
  semilogy(0:length(Nls{5})-1,Nls{5}, ...
           0:length(Nls{4})-1,Nls{4}, ...
           0:length(Nls{3})-1,Nls{3}, ...
           0:length(Nls{2})-1,Nls{2}, ...
           0:length(Nls{1})-1,Nls{1});
  xlabel('l'); ylabel('N_l'); %title(stitle)
  if option==4
    legend('\epsilon=0.0001','\epsilon=0.0002','\epsilon=0.0005','\epsilon=0.001','\epsilon=0.002',1)
  else
    legend('\epsilon=0.00005','\epsilon=0.0001','\epsilon=0.0002','\epsilon=0.0005','\epsilon=0.001',1)
  end
  axis([0 l 1e2 1e9])

  set(0,'DefaultAxesLineStyleOrder','-*|--*|-.*|:*')

  subplot(nvert,2,2*nvert)
  loglog(Eps,Eps.^2.*std_cost,Eps,Eps.^2.*mlmc_cost)
  xlabel('\epsilon'); ylabel('\epsilon^2 Cost'); %title(stitle)
  legend('Std MC','MLMC',1)
  if option==4
    axis([1e-4 2e-3 0.01 200])
    set(gca,'ytick',[0.01 0.1 1 10 100]);
  else
    axis([5e-5 1e-3 0.01 30])
  end

  if nvert==1
    if option==1
      print('-deps2c','mlmc2_gbm1b.eps')
    elseif option==2
      print('-deps2c','mlmc2_gbm2b.eps')
    elseif option==3
      print('-deps2c','mlmc2_gbm3b.eps')
    elseif option==4
      print('-deps2c','mlmc2_gbm4b.eps')
    elseif option==5
      print('-deps2c','mlmc2_gbm5b.eps')
    elseif option==6
      print('-deps2c','mlmc2_hestonb.eps')
    end

  else
    if option==1
      print('-deps2c','mlmc2_gbm1.eps')
    elseif option==2
      print('-deps2c','mlmc2_gbm2.eps')
    elseif option==3
      print('-deps2c','mlmc2_gbm3.eps')
    elseif option==4
      print('-deps2c','mlmc2_gbm4.eps')
    elseif option==5
      print('-deps2c','mlmc2_gbm5.eps')
    elseif option==6
      print('-deps2c','mlmc2_heston.eps')
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
B   = 0.85;

nf = M^l;
nc = nf/M;


hf = T/nf;
hc = T/nc;

if option==4
  nc = nc-1;
end

sums(1:6) = 0;

for N1 = 1:10000:N
  N2 = min(10000,N-N1+1);

%
% GBM model
%
  if option<6
    X0 = 1;

    Xf = X0*ones(1,N2);
    Xc = Xf;

    Af  = 0.5*hf*Xf;
    Ac  = 0.5*hc*Xc;

    Mf  = Xf;
    Mc  = Xc;
    Mc2 = Xc;

    Bf  = 1;
    Bc  = 1;
    Bc2 = 1;

    if l==0 & option~=4
      dWf = sqrt(hf)*randn(1,N2);
      Lf  = log(rand(1,N2));
      Xf0 = Xf;
      Xf  = Xf + r*Xf*hf + sig*Xf.*dWf + 0.5*sig^2*Xf.*(dWf.^2-hf);
      Af  = Af + 0.5*hf*Xf;
      vf  = sig*Xf0;
      Mf  = min(Mf,0.5*(Xf0+Xf-sqrt((Xf-Xf0).^2-2*hf*vf.^2.*Lf)));
      Bf  = Bf.*(1-exp(-2*max(0,(Xf0-B).*(Xf-B)./(hf*vf.^2))));
    else
      for n = 1:nc
        dWf = sqrt(hf)*randn(2,N2);
        Lf  = log(rand(3,N2));
        dIf = sqrt(hf/12)*hf*randn(2,N2);
        for m = 1:M
          Xf0 = Xf;
          Xf  = Xf + r*Xf*hf + sig*Xf.*dWf(m,:) + 0.5*sig^2*Xf.*(dWf(m,:).^2-hf);
          vf  = sig*Xf0;
          Af  = Af + hf*Xf + vf.*dIf(m,:);
          Mf  = min(Mf,0.5*(Xf0+Xf-sqrt((Xf-Xf0).^2-2*hf*vf.^2.*Lf(m,:))));
          Bf  = Bf.*(1-exp(-2*max(0,(Xf0-B).*(Xf-B)./(hf*vf.^2))));
        end

        dWc = dWf(1,:) + dWf(2,:);
        ddW = dWf(1,:) - dWf(2,:);

        Xc0 = Xc;
        Xc  = Xc + r*Xc*hc + sig*Xc.*dWc + 0.5*sig^2*Xc.*(dWc.^2-hc);

        vc  = sig*Xc0;
        Ac  = Ac + hc*Xc + vc.*(sum(dIf,1) + 0.25*hc*ddW);
        Xc1 = 0.5*(Xc0 + Xc + vc.*ddW);
        Mc  = min(Mc, 0.5*(Xc0+Xc1-sqrt((Xc1-Xc0).^2-2*hf*vc.^2.*Lf(1,:))));
        Mc  = min(Mc, 0.5*(Xc1+Xc -sqrt((Xc -Xc1).^2-2*hf*vc.^2.*Lf(2,:))));
        Mc2 = min(Mc2,0.5*(Xc0+Xc -sqrt((Xc -Xc0).^2-2*hc*vc.^2.*Lf(3,:))));
        Bc  = Bc .*(1-exp(-2*max(0,(Xc0-B).*(Xc1-B)./(hf*vc.^2))));
	Bc  = Bc .*(1-exp(-2*max(0,(Xc1-B).*(Xc -B)./(hf*vc.^2))));
	Bc2 = Bc2.*(1-exp(-2*max(0,(Xc0-B).*(Xc -B)./(hc*vc.^2))));
      end
      Af = Af - 0.5*hf*Xf;
      Ac = Ac - 0.5*hc*Xc;
    end

    if option==1
      Pf  = max(0,Xf-1);
      Pc  = max(0,Xc-1);
      Pc2 = Pc;
    elseif option==2
      Pf  = max(0,Af-1);
      Pc  = max(0,Ac-1);
      Pc2 = Pc;
    elseif option==3
      Pf  = Xf - Mf;
      Pc  = Xc - Mc;
      Pc2 = Xc - Mc2;
    elseif option==4
      if(l==0)
        Pf  = ncf((Xf+r*Xf*hf-1)./(sig*Xf*sqrt(hf)));
        Pc  = Pf;
        Pc2 = Pc;
      else
        dWf = sqrt(hf)*randn(1,N2);
        Xf  = Xf + r*Xf*hf + sig*Xf.*dWf + 0.5*sig^2*Xf.*(dWf.^2-hf);
        Pf  = ncf((Xf+r*Xf*hf-1)./(sig*Xf*sqrt(hf)));
        Pc  = ncf((Xc+r*Xc*hc+sig*Xc.*dWf-1)./(sig*Xc*sqrt(hf)));
        Pc2 = ncf((Xc+r*Xc*hc-1)./(sig*Xc*sqrt(hc)));
      end
    elseif option==5
      Pf  = Bf.*max(0,Xf-1);
      Pc  = Bc.*max(0,Xc-1);
      Pc2 = Bc2.*max(0,Xc-1);
    end

    dP  = exp(-r*T)*(Pf-Pc);
    dPc = exp(-r*T)*(Pc2-Pc);
    Pf  = exp(-r*T)*Pf;

%
% Heston model
%
  else
    r     =  0.05;
%    kappa =  5;
    kappa =  2;
    theta =  0.04;
    rho   = -0.5;
%    rho   = 0;
    xi    =  0.25;

    S0 = 1;
    v0 = theta;
    Sf = S0*ones(1,N2);
    vf = v0*ones(1,N2);
    Sc = Sf;
    vc = vf;
    dS = zeros(1,N2);

    if l==0
      dWf1 = sqrt(hf)*randn(1,N2);
      dWf2 = sqrt(hf)*randn(1,N2)*sqrt(1-rho^2) + rho*dWf1;
      Sf   = Sf.*(1 + r*hf + sqrt(max(0,vf)).*dWf1 ...
		  + 0.5*vf.*(dWf1.^2-hf) + 0.25*xi*(dWf1.*dWf2-rho*hf) );
      vf   = vf + (1-exp(-kappa*hf))*(theta-vf) ...
                 + exp(-kappa*hf)*xi*sqrt(max(0,vf)).*dWf2 ...
                  + exp(-kappa*hf)*0.25*xi^2*(dWf2.^2-hf);

    else
      Sf = [Sf Sf];
      vf = [vf vf];

      for n = 1:nc
        dWc1 = zeros(1,N2);
        dWc2 = zeros(1,N2);

        dWf1 = sqrt(hf)*randn(M,N2);
        dWf2 = sqrt(hf)*randn(M,N2)*sqrt(1-rho^2) + rho*dWf1;

        dWf1 = [dWf1 dWf1(M:-1:1,:)];
        dWf2 = [dWf2 dWf2(M:-1:1,:)];

        for m = 1:M
          Sf   = Sf.*(1 + r*hf + sqrt(max(0,vf)).*dWf1(m,:) ...
                      + 0.5*vf.*(dWf1(m,:).^2-hf) ...
                      + 0.25*xi*(dWf1(m,:).*dWf2(m,:)-rho*hf) );
          vf   = vf + (1-exp(-kappa*hf))*(theta-vf) ...
                     + exp(-kappa*hf)*xi*sqrt(max(0,vf)).*dWf2(m,:) ...
                      + exp(-kappa*hf)*0.25*xi^2*(dWf2(m,:).^2-hf);
        end

        dWc1  = dWf1(1,1:N2)+dWf1(2,1:N2);
        dWc2  = dWf2(1,1:N2)+dWf2(2,1:N2);

        Sc   = Sc.*(1 + r*hc + sqrt(max(0,vc)).*dWc1 ...
                    + 0.5*vc.*(dWc1.^2-hc) + 0.25*xi*(dWc1.*dWc2-rho*hc) );
        vc   = vc + (1-exp(-kappa*hc))*(theta-vc) ...
                   + exp(-kappa*hc)*xi*sqrt(max(0,vc)).*dWc2 ...
                    + exp(-kappa*hc)*0.25*xi^2*(dWc2.^2-hc);
      end
    end

    Pf = max(0,Sf-1);
    Pc = max(0,Sc-1);

    dPc = zeros(size(Pc));

    if l==0
      Pf = exp(-r*T)*Pf;
    else
      dP = exp(-r*T)*(0.5*(Pf(1:N2)+Pf(1+N2:end))-Pc(1:N2));
      Pf = exp(-r*T)*(0.5*(Pf(1:N2)+Pf(1+N2:end)));
    end
  end


  if l==0
    dP = Pf;
  end

  sums(1) = sums(1) + sum(dP);
  sums(2) = sums(2) + sum(dP.^2);
  sums(3) = sums(3) + sum(Pf);
  sums(4) = sums(4) + sum(Pf.^2);
  sums(5) = sums(5) + sum(dPc);
  sums(6) = sums(6) + sum(dPc.^2);
end

