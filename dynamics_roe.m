function xdot = dynamics_roe(t, x,n)

A = zeros(6,6); A(2,1) = -3/2*n; 
xdot = A*x; 