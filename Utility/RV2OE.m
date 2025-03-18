
function [a,e,INC,AOP,RAAN,MA,TA] = RV2OE(R,V,mu)
%--------------------------------------------------------------------------
%   Converts R and V to classical orbital elements
%--------------------------------------------------------------------------
%   Form:
%   [a,e,INC,AOP,RAAN,MA,TA] = RV2OE( R, V, mu)
%--------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   R               (1,3) Position vector
%   V               (1,3) Velocity vector
%   mu                    Gravitational parameter
%
%   -------
%   Outputs
%   -------
%   X               (1,6) Elements vector [a,e,INC,AOP,RAAN,MA,TA]
%
% Magnitude of vectors (lowercase)
r = norm(R);    v = norm(V);
% Calculate h and e
H = cross(R,V);     % angular momentum
h = norm(H);
E = -R/r - 1/mu*cross(H,V);   % eccentricity
e = norm(E);
% Perifocal frame
u1 = E/e;
u3 = H/h;
u2 = cross(u3,u1);
% Inertial frame
i = [1 0 0];
j = [0 1 0];
k = [0 0 1];
uln = cross(k,u3);    % node line
% Angles from perifocal frame
INC = acos(dot(u3,k));    % inclination
RAAN = atan2(dot(uln,j),dot(uln,i));    % ascending node
if RAAN < 0
    RAAN = RAAN + 2*pi;
end
AOP = atan2(dot(-u2,uln),dot(uln,u1));  % argument of perigee
if AOP < 0
    AOP = AOP + 2*pi;
end
TA = atan2(dot(R,u2),dot(R,u1));    % true anomaly
if TA < 0
    TA = TA + 2*pi;
end
u = 2*atan(sqrt((1-e)/(1+e))*tan(TA/2));    % eccentric anomaly
MA = u - e*sin(u);  % mean anomaly (Kepler's eq.)
% Add semimajor axis from energy equation
a = r/(2 - v^2*r/mu);
end
