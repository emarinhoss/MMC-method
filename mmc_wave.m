% 1D-wave equation
function sums = mmc_wave(M,l,N,T,R,nr)

mf = 2^l*nr;
mc = 2^(l-1)*nr;

nf = M^l;
nc = nf/M;

% grid
rf = linspace(-R,R,mf+1); df = abs(rf(2)-rf(1));
rc = linspace(-R,R,mc+1); dc = abs(rc(2)-rc(1));
% initial condition
wf = exp(-rf.^2);
wc = exp(-rc.^2);

cfl = 1.0;
rr=linspace(-R,R,nr+1);
sums = zeros(4,length(rr));
for N1 = 1:1e4:N
    N2 = min(1e4,N-N1+1);
    
    c = -1 + 0.1*randn(1,N2);
    dtf = abs(cfl*df./c);
    dtc = abs(cfl*dc./c);
    
    solf = zeros(N2,length(rf)+1);    
    solc = zeros(N2,length(rc)+1);
    
    %%%
    if l==0
        for m=1:N2
            solf(m,1:end-1)=wf;        % IC
            solf(m,end)=solf(m,end-1); % BC-up
            t = 0;
            while (t <= T)
                solf(m,1:end-1) = solf(m,1:end-1) - ...
                    (dtf(m)/df)*c(m)*(solf(m,2:end) - solf(m,1:end-1));
                solf(m,end) = solf(m,end-1);
                t = t + dtf(m);
            end
        end
    else
        for m=1:N2
            solf(m,1:end-1)=wf;        % IC
            solf(m,end)=solf(m,end-1); % BC-up
            t = 0;
            while (t <= T)
                solf(m,1:end-1) = solf(m,1:end-1) - ...
                    (dtf(m)/df)*c(m)*(solf(m,2:end) - solf(m,1:end-1));
                solf(m,end) = solf(m,end-1);
                t = t + dtf(m);
            end
            
            solc(m,1:end-1)=wc;       % IC
            solc(m,end)=solc(m,end-1);% BC-up
            t = 0;
            while (t <= T)
                solc(m,1:end-1) = solc(m,1:end-1) - ...
                    (dtc(m)/df)*c(m)*(solc(m,2:end) - solc(m,1:end-1));
                solc(m,end) = solc(m,end-1);
                t = t + dtc(m);
            end
        end
    end
    
    Pf = zeros (1,length(rr)); Pc = Pf; Pfs = Pf; Pcs = Pf;
    for m=1:N2
        f = interp1(-R:df:R,solf(m,1:end-1),rr);
        c = interp1(-R:dc:R,solc(m,1:end-1),rr);
        
        Pf = Pf + f; Pfs = Pfs + f.*f;
        Pc = Pc + c; Pcs = Pcs + c.*c;
    end
    
    sums(1,:) = sums(1,:) + (Pf - Pc);
    sums(2,:) = sums(2,:) + (Pfs - Pcs);
    sums(3,:) = sums(3,:) + (Pf);
    sums(4,:) = sums(4,:) + (Pfs);
end