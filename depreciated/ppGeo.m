
%param.TU = 86400;  % s
mu = 3.986005e14;%m^2/s^3
Re = 6378137; % m
J2 =  1.08262668e-3;
N = 25;
k = 3*J2*Re^2/(2*mu^3); %
AU =  149597870691; %m
Rs = 695800e3; %m
g0 = 9.80665; %m/s
m_servicer_wet = 800;


param.LU = Re;
param.VU = sqrt(mu/Re);
param.TU = param.LU/param.VU;
param.MU = m_servicer_wet;


param.mu= mu/param.LU^3*param.TU^2; 
param.Re = Re/param.LU;
param.J2 = J2; 
param.N = N;
param.k = k/param.TU^6*param.LU^7;
param.AU = AU/param.LU;
param.Rs = Rs/ param.LU;
param.g0 = g0/param.LU*param.TU*param.TU; %m/s

%% Problem dynamics
param.eclipses = false; % consider eclipses.
param.drag = false;
param.dutyRatio =1; % consider dutyRatio
param.Isp = 1300/param.TU; %s
param.g0 = 9.80665/param.LU*param.TU*param.TU; %m/s

param.IspImp = 300/param.TU;
param.TImp = 1/param.MU/param.LU*param.TU^2;

param.a_pchangeLimit = 0;
param.Tmax =60e-3/param.MU/param.LU*param.TU^2;
param.Cd = 2.2;
param.Area= 5/param.LU/param.LU;
param.dragfactor = 1; 
param.m_servicer_wet = m_servicer_wet/param.MU; 
param.target_altitude = 350e3/param.LU;

param.m0  = param.m_servicer_wet;

datevec = '19-June-2023';
param.t0 = juliandate(datevec,'dd-mmm-yyyy'); %+ 1000;
param.t0JD = param.t0; 
eci_t = (param.t0 - juliandate('01-January-2000','dd-mmm-yyyy'))*24*60*60;

[rr, vv, debris(1,:), m_debris(1)] = getPosition(param.t0, ...
       'H2AF15',3.986005e14);
   eci_t = (param.t0 - juliandate('01-January-2000','dd-mmm-yyyy'))*24*60*60;

eci_pos = [rr;vv]/1e3;

a0 = 38500e3/param.LU + param.Re;
inc0 = debris(1,3) + pi/180;
RAAN0_t0 =debris(1,4);
omega0 =  debris(1,5);
e0 =  debris(1,2);
v0 = sqrt(param.mu/a0);
