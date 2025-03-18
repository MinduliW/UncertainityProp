
function [ nu ] = eclipse( x_sun, x_sat, param)
%Eclipse function: calculates the degree of eclipse due to Earth
%From Montebruck and Gill algorithm
%INPUT:
%       - x_sun: 3D vector of the Sun position and velocity in Inertial ref. frame
%       - x_sat: 3D vector of the satellite position and velocity in Inertial ref. frame
% OUTPUT:
%        -nu: eclipse value (0 in eclipse, 1 in sunligth, (0,1) in penumbra)


r_sat = norm(x_sat);

a = asin(param.Rs/ norm(x_sun - x_sat));
b = asin(param.Re/r_sat);
c = acos(dot((x_sun-x_sat)/norm(x_sun - x_sat),-x_sat/r_sat));

if a+b<= c
    nu = 1;
elseif c <= abs(b-a)
    nu = 0;
else
    x = (c^2 + a^2 - b^2 ) / (2*c);
    y = sqrt(a^2 - x^2);
    
    A = a^2 * acos(x/a) + b^2 * acos((c-x)/b)- c*y;
    
    nu = 1 - A/(pi*a^2);
end


