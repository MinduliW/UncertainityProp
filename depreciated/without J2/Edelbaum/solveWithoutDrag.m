function [res,ad,RAAN,wtg]  = solveWithoutDrag...
    (waitTime,Vd,Id, RAAN0,massAfterL1, timeAfterL1, param)

% output : residual on final position acheived. 
% input : wait time 

waitTime = (waitTime^2*param.TU);
% calculate the drift orbit altitude. 
adrift = param.mu/Vd^2;
% calculate the drift orbit period. 
Td = 2*pi*sqrt(adrift^3/param.mu);

% divide waittime to orbital periods.
wtg = 0:Td:waitTime; 


ad = adrift*ones(size(wtg)); %m

        
Vd = sqrt(param.mu./ad);
        
RAAN =  RAAN0-param.k*Vd(1).^7*cos(Id)*(wtg);
        

% calculate the effect of the second leg. 

param.m0 = massAfterL1; 

param.timeprev = timeAfterL1+waitTime; 

af = param.mu/param.xf(1)^2;

[t2,~, Omega_tf] =...
    kluver(ad(end),af, Id,param.xf(2) ,RAAN(end),param);

% target RAAN
TOF = param.timeprev + t2(end); 
targetRAAN = param.xf(3) -param.k*param.xf(1)^7*cos(param.xf(2))*(TOF);



res = sin(0.5*(targetRAAN - Omega_tf));


end
