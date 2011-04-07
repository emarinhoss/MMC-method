function dhdt = tankfill(t,h,dummy,gamma,beta)
% RHS funtion for tank-fill problem
% dh/dt = gamma*alpha(t) - beta*sqrt(h), h(0) = ho

alpha = sin(t);

dhdt = 10 + gamma*alpha - beta*sqrt(h);