%
% function V = heston_call(r,kappa,theta,omega,rho,S,K,V0,T)
%
% analytic solution for European call option in Heston model, using
% algorithm in "Not-so-complex logarithms in the Heston model"
% by Christian Kahl and Peter Jaekel
% http://www.math.uni-wuppertal.de/~kahl/publications/NotSoComplexLogarithmsInTheHestonModel.pdf
% http://www.btinternet.com/~pjaeckel/NotSoComplexLogarithmsInTheHestonModel.pdf
%
% r     - interest rate
% kappa - mean reversion rate
% theta - mean reversion volatility
% omega - vol-of-vol coefficient
% rho   - correlation factor
% S     - asset value(s)
% K     - strike price
% V0    - initial volatility
% T     - time interval
% V     - option value(s)

function V = heston_call(r,kappa,theta,omega,rho,S,K,V0,T)

if nargin ~= 9
  error('wrong number of arguments');
end

V = quadl(@(x)integrand(x,r,kappa,theta,omega,rho,S,K,V0,T),0,1,1e-12);
V = exp(-r*T)*V;

%--------------------------------------------------------

function y = integrand(x,r,kappa,theta,omega,rho,S,K,V0,T)

x    = max(1e-20,min(x,1-1e-10));
Cinf = sqrt(1-rho^2)*(V0+kappa*theta*T)/omega;
u    = - log(x)/Cinf;
F    = S*exp(r*T);

for ipass = 1:2
  um  = u - i*(2-ipass);
  d   = sqrt( (rho*omega*um*i - kappa).^2 + omega^2*(um*i + um.^2));
  c   = (kappa-rho*omega*um*i+d)./(kappa-rho*omega*um*i-d);
  tc  = angle(c);
  GD  = c-1;
  m   = floor( (tc+pi)/(2*pi) );
  GN  = c.*exp(i*imag(d)*T) - exp(-real(d)*T);
  n   = floor( (tc+imag(d)*T+pi)/(2*pi) );
  lnG = real(d)*T + log(abs(GN)./abs(GD)) + i*(angle(GN)-angle(GD)+2*pi*(n-m));
  D   = ((kappa-rho*omega*um*i+d)/omega^2) .* ((exp(d*T)-1)./(c.*exp(d*T)-1));
  C   = ((kappa*theta)/omega^2) * ( (kappa-rho*omega*um*i+d)*T - 2*lnG );
  phi = exp(C+D*V0+i*um*log(F));
  f(ipass,:) = real( exp(-i*u*log(K)).*phi ./ (i*u) );
end

y = 0.5*(F-K) + (f(1,:) - K*f(2,:))./(x*pi*Cinf);

