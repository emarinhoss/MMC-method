function vals = MMC_tank(M,l,N)
% integration interval
T = 1;

nf = M^l;
nc = nf/M;

hf = T/nf;
hc = T/nc;

vals(1:2) = 0;

for N1 = 1:10000:N
    %
    N2 = min(10000,N-N1+1);
    
    % Initial height of liquid in the tank
    h0 = 1;
    
    hl = h0*ones(1,N2);
    hlm1 = hl;
    
    t =0;
    tt=0;
    beta = randn(1,N2);
    if l==0
      hl  = hl + hf*(10+4*sin(t)-beta.*sqrt(hl));
    else
      for n = 1:nc
        for m = 1:M
          hl  = hl + hf*(10+4*sin(t)-beta.*sqrt(hl));
          t=t+hf;
        end
        hlm1 = hlm1 + hc*(10+4*sin(tt)-beta.*sqrt(hlm1));
        tt=tt+hc;
      end
    end
end

vals(1) = sum(hl-hlm1);
vals(2) = sum((hl-hlm1).*(hl-hlm1));
