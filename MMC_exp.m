function vals = MMC_exp(M,l,N,T,y0)

nf = M^l;
nc = nf/M;

hf = T/nf;
hc = T/nc;

vals(1:2) = 0;

for N1 = 1:10000:N
    
    N2 = min(10000,N-N1+1);
        
    hl = y0*ones(1,N2);
    hlm1 = hl;
    
    
    if l==0
      k = randn(1,N2);
      hl  = hl - hf*(k.*hl);
    else
      for n = 1:nc
        for m = 1:M
          k = randn(1,N2);
          hl  = hl - hf*(k.*hl);
        end
        k = randn(1,N2);
        hlm1 = hlm1 - hc*(k.*hlm1);
      end
    end
end

vals(1) = sum(hl-hlm1);
vals(2) = sum((hl-hlm1).*(hl-hlm1));
