%
% function V = digital_call(r,sigma,T,S,opt)
%
% Black-Scholes digital call option solution
% as defined on page 82 of 
% The Mathematics of Financial Derivatives
% by Wilmott, Howison and Dewynne
%
% r     - interest rate
% sigma - volatility
% T     - time interval
% S     - asset value(s) / strike price
% opt   - 'value', 'delta', 'gamma' or 'vega'
% V     - option value(s)
%

function V = digital_call(r,sigma,T,S,opt)

S  = max(1e-20,S);

d2 = ( log(S) + (r-0.5*sigma^2)*T ) / (sigma*sqrt(T));

switch lower(opt)
  case 'value'
    V = exp(-r*T)*N(d2);

  case 'delta'
    V = exp(-r*T-d2.^2/2) ./ (S*sigma*sqrt(2*pi*T));

  case 'gamma'
    V = exp(-r*T-d2.^2/2) .* ( -d2./(sqrt(2*pi)*S.^2*sigma^2*T) ...
                               - 1./(sqrt(2*pi*T)*S.^2*sigma) );
  case 'vega'
    V = exp(-r*T-d2.^2/2) .* (-sqrt(T)-d2/sigma) ./ sqrt(2*pi);

  otherwise
    error('invalid value for opt -- must be value, delta or gamma')
end





function cum = N(x)

%cum = 0.5*(1+erf(x/sqrt(2)));

xr = real(x);
xi = imag(x);

if abs(xi)>1e-10
  error 'imag(x) too large in N(x)'
end

cum = 0.5*(1+erf(xr/sqrt(2))) ...
    + 1i*xi.*exp(-0.5*xr.^2)/sqrt(2*pi);
