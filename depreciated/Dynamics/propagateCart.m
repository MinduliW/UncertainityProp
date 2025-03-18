function xdot = propagateCart(t,x,param)
 
% FUNCTION NAME:
%   propagateCart
%
% DESCRIPTION:
%   Propagates dynamics, USed for testing only.
%
% INPUT:
%   t - (double) time
%   x - (double []) x in Cartesian coordinates 
%   param - (struct) Problem parameters 
%
% OUTPUT:
%   xdot - (double []) time derivative of state x. 
%


r = sqrt(x(1)^2+ x(2)^2 + x(3)^2);

xdot(1) = x(4); 
xdot(2) = x(5); 
xdot(3) = x(6); 
xdot(4) = -param.mu*x(1)/r^3; 
xdot(5) = -param.mu*x(2)/r^3; 
xdot(6) = -param.mu*x(3)/r^3; 


xdot(4) =xdot(4) -param.mu*x(1)/r^3*1.5*param.J2*(param.Re/r)^2*(1- 5*x(3)^2/r^2);
xdot(5) =xdot(5)-param.mu*x(2)/r^3*1.5*param.J2*(param.Re/r)^2*(1- 5*x(3)^2/r^2);
xdot(6) =xdot(6) -param.mu*x(3)/r^3*1.5*param.J2*(param.Re/r)^2*(3- 5*x(3)^2/r^2);

% introduce drag
rho = Density_HP((r - param.Re)*param.LU)*param.LU*param.LU*param.LU/param.MU;


v=  sqrt(x(4)^2 + x(5)^2 + x(6)^2);

vhat = x(4:6)/v; 

if param.drag == true
    drag = 0.5*rho*param.Cd *param.Area*v^2/param.m0*(-vhat);
else
    drag = [0;0;0];
end

xdot(4) = xdot(4)+ drag(1);
xdot(5) = xdot(5)+ drag(2);
xdot(6) = xdot(6)+ drag(3);


xdot = xdot';
end
