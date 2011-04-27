function sums = MMC_exp(M,l,N,T,y0,k1,k2)

nf = M^l;
nc = nf/M;

hf = T/nf;
hc = T/nc;

sums(1:4) = 0;
for N1 = 1:1e5:N
    N2 = min(1e5,N-N1+1);
    dWf = k1 + k2*randn(1,N2);

    Xf = y0*ones(1,N2);
    Xc = Xf;

    if l==0
      Xf  = Xf + hf * (Xf.*dWf);
    else
      for n = 1:nc
        for m = 1:M
          Xf  = Xf + hf * (Xf.*dWf);
        end
        Xc = Xc + hc * (Xc.*dWf);
      end
    end

  if l==0
    Xc=0;
  end

  sums(1) = sums(1) + sum(Xf-Xc);
  sums(2) = sums(2) + sum((Xf-Xc).*(Xf-Xc));
  sums(3) = sums(3) + sum(Xf);
  sums(4) = sums(4) + sum(Xf.*Xf);
end