function gamma_L = Smoothed_eclipse(cs,ct, x_sun, x_sat, param)

% Radius of Earth and Sun.
Rs = param.Rs; 
Re = param.Re; 

a_SR = asin(Rs / norm(x_sun)); 
a_BR = asin(Re/norm(x_sat)); 

a_D = acos(dot(x_sat,x_sun)/norm(x_sat)/norm(x_sun));

gamma_L = 1/(1+ exp(-cs*(a_D -ct*(a_SR + a_BR))));

end