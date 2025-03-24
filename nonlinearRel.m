function xdot = nonlinearRel(t,x,param)
 
a = param.a;
mu = param.mu; 


ax = 0;
ay = 0;
az = 0; 



n = sqrt(mu/a^3); 

xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);

r = sqrt((x(1)+a)^2 + x(2)^2 + x(3)^2); 

xdot(4) = ax + 2*n*x(5) +n^2*x(1) + n^2*a - mu*(x(1)+a)/(r^3);
xdot(5) = ay + n^2*x(2) - 2*n*x(4) -mu*x(2)/(r^3); 
xdot(6) = az + -mu*x(3)/((r)^3);

xdot = xdot';

