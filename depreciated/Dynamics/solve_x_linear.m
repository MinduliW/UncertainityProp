function xdot = solve_x_linear(t,xvalues, param)
% FUNCTION NAME:
%   solve_x_linear
%
% DESCRIPTION:
%   Describes the linearised energy optimal dynamics 
%
% INPUT:
%   t - (double) time
%   x - (double []) x in MEE coordinates 
%   param - (struct) Problem parameters 
%
% OUTPUT:
%   xdot - (double []) time derivative of state x. 
%

x = param.x0 + (param.xf - param.x0)/(param.tf - param.t0)*(t - param.t0);

% Generate new A B C D matricies.
[~,B,dAdx, D, ~] = getAandBMatrices(x',param);

Blambda = B'*xvalues(8:13);

if param.eclipse == true
    %Shadow function computation
   v = 1- getEclipse(t, xvalues, param);

    
else
    v = 1; 
    
end

xdot(1:6,1) =  dAdx*xvalues(1:6) + D - v*param.Tmax/xvalues(7)*B*Blambda;

xdot(7,1) = -v*param.Tmax*norm(Blambda)/(param.Isp*param.g0);

xdot(8:13,1) =  -(dAdx)'*xvalues(8:13);



if param.addJ2 == true
    % Get J2 terms
    [deltaJ2,~] = J2acceleration(xvalues, param);
    deltaJ2 = deltaJ2';
    xdot(1:6,1) =  xdot(1:6,1) + B*deltaJ2;
end




end