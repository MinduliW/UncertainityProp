
%param.TU = 86400;  % s
mu = 3.986005e14;%m^2/s^3
Re = 6378.137e3; % m
J2 =  1.08262668e-3;

g0 = 9.80665; %m/s
m_servicer_wet = 800;


param.LU = Re;
param.VU = sqrt(mu/Re);
param.TU = param.LU/param.VU;
param.MU = m_servicer_wet;



param.mu= mu/param.LU^3*param.TU^2; 
param.Re = Re/param.LU;
param.J2 = J2; 


a0 = 38500e3/param.LU + param.Re;
