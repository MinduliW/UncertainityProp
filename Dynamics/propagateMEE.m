function xdot = propagateMEE(t,x,param)
 
% FUNCTION NAME:
%   propagateMEE
%
% DESCRIPTION:
%   Propagates dynamics with the control obtained from the initial guess.
%   USed for testing only.
%
% INPUT:
%   t - (double) time
%   x - (double []) x in MEE coordinates 
%   param - (struct) Problem parameters 
%
% OUTPUT:
%   xdot - (double []) time derivative of state x. 
%

% x = param.x0 + (param.xf - param.x0)/(param.tf - param.t0)*(t - param.t0);
% 
% % Generate new A B C D matricies.
% [~,B,dAdx, D, ~] = getAandBMatrices(x',param);
% 
% 
% ux = interp1(param.time, param.control(:,1), t); 
% uy = interp1(param.time, param.control(:,2), t); 
% uz = interp1(param.time, param.control(:,3), t); 
% 
% u = [ux;uy;uz];
% 
% xdot(1:6,1) =  dAdx*xvalues(1:6) + D +B*u;
% 
% xdot(7,1) = 0;


%Generate new A B C D matricies.
[A,B] = getAandBMatrices(x,param);


ux = interp1(param.time, param.control(:,1), t);
uy = interp1(param.time, param.control(:,2), t);
uz = interp1(param.time, param.control(:,3), t);

u = [ux;uy;uz];

xdot(1:6,1) =  A + B*u;

xdot(7,1) = -param.Tmax*norm(u)/(param.Isp*param.g0);

end


