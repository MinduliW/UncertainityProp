function X = RV2GEq(R,V,mu,J2, Re)
%--------------------------------------------------------------------------
%   Converts R and V to Equinoctial elements
%--------------------------------------------------------------------------
%   Form:
%   eq = RV2GEq( R, V, mu ,U)
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   R               (1,3) Position vector
%   V               (1,3) Velocity vector
%   mu                    Gravitational parameter
%   U               (1,1) Perurbation potential (i.e. negative disturbing
%                         function). For U=0 classical "alternate" 
%                         equinoctial elements are treated.
%
%   -------
%   Outputs
%   -------
%   X               (1,6) Elements vector [nu,p1,p2,varL,q1,Q2]
%
%--------------------------------------------------------------------------
%   References:  Baù, G., Hernando-Ayuso, J. and Bombardelli, C., 2021.
%   "A generalization of the equinoctial orbital elements" 
%   Celestial Mechanics and Dynamical Astronomy, 133(11), pp.1-29.
%--------------------------------------------------------------------------


%% Calculate GEqOE from R,V
[~,~,INC,~,RAAN,~,~] = RV2OE(R,V,mu);  % obtain RAAN and INC
% State parameters
r = norm(R);
v = norm(V);
h = norm(cross(R,V));           % angular momentum module
r_dot = dot(R,V)/r;             % radial velocity
% Calculate energy
epsK = 1/2*v^2 - mu/r;          % Keplerian energy

phi = acos(R(3)/r); 
U = J2/2*mu/r*(Re/r)^2*(3*cos(phi)^2-1);

eps = epsK + U;                 % total energy
nu = 1/mu*(-2*eps)^(3/2);       % (1)
% Calculate q in order to obtain equinoctial frame
q1 = tan(INC/2)*sin(RAAN);  % (5)
q2 = tan(INC/2)*cos(RAAN);  % (6)
% Unit vectors of equinoctial frame in inertial ref.
Kq = 1/(1+q1^2+q2^2);
ex = Kq*[1-q1^2+q2^2    ,2*q1*q2        ,-2*q1];
ey = Kq*[2*q1*q2        ,1+q1^2-q2^2    ,2*q2];
% Radial vector in inertial ref.
er = R/r;
% Obtain true longitude
cL = dot(er,ex);
sL = dot(er,ey);
L = atan2(sL,cL);
if L < 0
    L = L + 2*pi;
end
% L = AOP + RAAN + TA;
% Calculate rho,c and a
Ueff = h^2/2/r^2 + U;   % effective potential
c = sqrt(2*r^2*Ueff);
rho = c^2/mu;
a = -mu/2/eps;
% g vector
p1 = (rho/r-1)*sin(L) - c*r_dot/mu*cos(L);  % (2)
p2 = (rho/r-1)*cos(L) + c*r_dot/mu*sin(L);  % (3)
% Obtain K for varL
w = sqrt(mu/a);
sK = (mu+c*w-r*r_dot^2)*sin(L) - r_dot*(c+w*r)*cos(L);
cK = (mu+c*w-r*r_dot^2)*cos(L) + r_dot*(c+w*r)*sin(L);
K = atan2(sK,cK);
if K < 0
    K = K + 2*pi;
end
% Obtain generalized mean longitud varL
varL = K + 1/(mu+c*w)*(cK*p1-sK*p2);  % (4)
% GEqOE array
X = [nu p1 p2 varL q1 q2];
end

