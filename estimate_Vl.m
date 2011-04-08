function values = estimate_Vl(hl,N,tf)

% height
h = ones(N,length(hl),iter);

gamma = randn(1,N);
beta  = randn(1,N);

for k=1:N
    for i=1:length(hl)
        t = 0;
        for j=1:iter
            t=t+hl;
            h(k,i,j) = h(k,i,j) + hl(i) .* (10+gamma(k).*sin(t)-...
                beta(k).*sqrt(h(k,i,j)));
        end
    end
end
