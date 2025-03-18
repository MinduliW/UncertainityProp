function [R,V] = GEq2RV(X,tol,mu,J2, Re)
%--------------------------------------------------------------------------
%   Converts generalized equinoctial elements to r and v 
%   for an elliptic orbit.
%--------------------------------------------------------------------------
%   Form:
%   [r, v] =  GEq2RV( X, tol, mu, U )
%--------------------------------------------------------------------------
%
%   ------ 
%   Inputs
%   ------
%   X               (1,6) Elements vector [nu,p1,p2,varL,q1,Q2]
%   tol             (1,1) Tolerance for Kepler's equation solver
%   mu              (1,1) Gravitational parameter
%   U               (1,1) Perurbation potential (i.e. negative disturbing
%                         function). For U=0 classical "alternate" 
%                         equinoctial elements are treated.
%
%   -------
%   Outputs
%   -------
%   R               (1,3) position vector
%   V               (1,3) velocity vector
%
%--------------------------------------------------------------------------
%   References:  Baù, G., Hernando-Ayuso, J. and Bombardelli, C., 2021.
%   "A generalization of the equinoctial orbital elements" 
%   Celestial Mechanics and Dynamical Astronomy, 133(11), pp.1-29.
%--------------------------------------------------------------------------

nu = X(1); p1 = X(2); p2 = X(3); varL = X(4); q1 = X(5); q2 = X(6);
%% Auxiliary variables
a = (mu/nu^2)^(1/3);                    % semi-major axis
c = (mu^2/nu)^(1/3)*sqrt(1-p1^2-p2^2);  % generalized angular momentum
alpha = 1/(1 + sqrt(1-p1^2-p2^2));
%% Obtain true longitude L
% Solve generalized Kepler eq. for K
fun = @(K)gen_kepler(K,p1,p2,varL);
opts = optimset('TolX',tol);
K = fzero(fun,varL,opts);
% r and r_dot once known K
r = a*(1-p1*sin(K)-p2*cos(K));              % radius [LU]



r_dot = sqrt(mu*a)/r*(p2*sin(K)-p1*cos(K)); % radius time der [LU/TU]
% Calculate L
sL = a/r*(alpha*p1*p2*cos(K)+(1-alpha*p2^2)*sin(K)-p1);     % sin
cL = a/r*(alpha*p1*p2*sin(K)+(1-alpha*p1^2)*cos(K)-p2);     % cos
L = atan2(sL,cL);
if L < 0
    L = L + 2*pi;
end
%% Frame rotation and position vector
% inertial to equinoctial
Kq = 1/(1+q1^2+q2^2);
ex = Kq*[1-q1^2+q2^2    ,2*q1*q2        ,-2*q1];
ey = Kq*[2*q1*q2        ,1+q1^2-q2^2    ,2*q2];
% equinoccial to local
er = ex*cos(L) + ey*sin(L);
ef = ey*cos(L) - ex*sin(L);  
% Position components in equinoctial frame
R = r*er;    % position vector in inertial frame

phi = acos(R(3)/r); 
U = J2/2*mu/r*(Re/r)^2*(3*cos(phi)^2-1);


%% Velocity vector
% Variables that required perturbations
h = sqrt(c^2 - 2*r^2*U);    % angular moment. module [km^2/s]
V = r_dot*er + h/r*ef;      % velocity vector in inertial frame
end

