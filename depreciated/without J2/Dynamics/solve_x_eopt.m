function xdot = solve_x_eopt(t,x, param)

% FUNCTION NAME:
%   solve_x_eopt
%
% DESCRIPTION:
%   Describes the nonlinear energy optimal dynamics 
%
% INPUT:
%   t - (double) time
%   x - (double []) x in MEE coordinates 
%   param - (struct) Problem parameters 
%
% OUTPUT:
%   xdot - (double []) time derivative of state x. 
%
% Generate new A B C D matricies.
[A,B,dAdx] = getAandBMatrices(x,param);

Blambda = B'*x(8:13);

xdot(1:6,1) =  A - param.Tmax/x(7)*B*Blambda;

xdot(7,1) = 0; %-param.Tmax*norm(Blambda)/(param.Isp*param.g0);

xdot(8:13,1) =  -(dAdx)'*x(8:13);


end