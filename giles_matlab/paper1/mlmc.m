% function [P, Nl] = mlmc(M,eps,mlmc_l)
%
% multi-level Monte Carlo path estimation
%  P      = value
%  Nl     = number of samples at each level
%  M      = timestep refinement factor
%  eps    = accuracy (rms error)
%  mlmc_l = function for level l estimator 
%
%  mlmc_l(M,l,N)
%       l = level
%       N = number of paths

function [P, Nl] = mlmc(M,eps,mlmc_l,extrap)

L   = -1;
N   = 10000;
converged = 0;

while ~converged
%
% initial variance estimate
%
  L = L+1;
  sums = feval(mlmc_l,M,L,N);
  suml(1,L+1) = N;
  suml(2,L+1) = sums(1);
  suml(3,L+1) = sums(2);

%
% optimal sample sizes
%
  Vl = suml(3,:)./suml(1,:) - (suml(2,:)./suml(1,:)).^2;
%  Nl = ceil( 2 * Nl * sum(Vl./Nl) / eps^2);
  Nl = ceil( 2 * sqrt(Vl./(M.^(0:L))) * sum(sqrt(Vl.*(M.^(0:L)))) / eps^2);
%  disp(sprintf(' %d ',Nl))

%
% update sample sums
%
  for l=0:L
    dNl = Nl(l+1)-suml(1,l+1);
    if dNl>0
      sums = feval(mlmc_l,M,l,dNl);
      suml(1,l+1) = suml(1,l+1) + dNl;
      suml(2,l+1) = suml(2,l+1) + sums(1);
      suml(3,l+1) = suml(3,l+1) + sums(2);
    end
  end

%
% test for convergence
%

  if extrap
    range = 0;
    if L>1 & M^L>=16
      con = M.^range.*( suml(2,L+1+range)./suml(1,L+1+range) ...
                - (1/M)*suml(2,L  +range)./suml(1,L  +range) );
      converged = (max(abs(con)) < (M^2-1)*eps/sqrt(2)) | (M^L>1024) ;
    end
  else
    range = -1:0;
    if L>1 & M^L>=16
      con = M.^range.*suml(2,L+1+range)./suml(1,L+1+range);
%      converged = (max(abs(con)) < (M-1)*eps/2) | (M^L>1024) ;
      converged = (max(abs(con)) < (M-1)*eps/sqrt(2)) | (M^L>1024) ;
    end
  end
end

%
% evaluate multi-timestep estimator
%
  P = sum( suml(2,1:L+1)./suml(1,1:L+1) );
  if extrap
    P = P + ( suml(2,L+1)./suml(1,L+1) ) / (M-1);
  end
end
