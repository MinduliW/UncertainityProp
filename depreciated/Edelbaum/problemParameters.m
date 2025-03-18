constants; 

%% Data
param.eclipses = false; % consider eclipses.
param.drag = false;
param.dutyRatio = 1; % consider dutyRatio
param.Isp = 1300; %s
param.g0 = 9.80665; %m/s

param.a_pchangeLimit = 0;
param.T =60e-3;
param.Cd = 2.2;
param.Area= 8;
param.dragfactor = 1; 
param.m_servicer_wet = 800; 
param.target_altitude = 350e3;

param.m0  = param.m_servicer_wet;

datevec = '19-June-2024';
param.t0 = juliandate(datevec,'dd-mmm-yyyy'); %+ 1000;

debrisID = {'H2AF15','H2AF24'};
for i = 1:length(debrisID)
    [~, ~, debris(i,:), ~, m_debris(i)] = getPosition(param.t0 , ...
        debrisID(i),param.mu, param.J2, param.Re);
    
end


a0 = param.target_altitude + param.Re; %+90e3; %+150e3;
inc0 = debris(1,3);
RAAN0_t0 =debris(1,4);
v0 = sqrt(param.mu/a0);

param.x0 = [v0,  inc0, RAAN0_t0]; 


af = debris(2,1);
incf = debris(2,3)-pi/2/180;
RAANf_t0 =  debris(2,4);
vf = sqrt(param.mu/af);

param.xf = [vf , incf, RAANf_t0]; 

