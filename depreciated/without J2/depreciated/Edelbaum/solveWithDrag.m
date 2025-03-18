function [res,ad,RAAN,wtg,decaystatus]  = solveWithDrag(waitTime,Vd,Id, RAAN0,massAfterL1, timeAfterL1, param)

% output : residual on final position acheived.
% input : wait time

waitTime = (waitTime^2*param.TU);

% calculate the drift orbit altitude.
adrift = param.mu/Vd^2;

% density evaluated at perigee
rho_hp = Density_HP(adrift - param.Re);

% area to mass ratio
delta = param.Area/param.m0;

t = 0;i = 1;
wtg = [];
wtg(1) = t; 
ad(1) = adrift;
RAAN(1) = RAAN0;

while t < waitTime
    
    % from Extension of the King-Hele orbit contraction method for accurate, semi-analytical propagation of non-circular orbit
    delta_a = -2*pi*delta*ad(i)^2*rho_hp;
    

   ad(i+1) = ad(i)+ delta_a;
    

    Vdrift = 0.5*(sqrt(param.mu/ad(i+1))+ sqrt(param.mu/ad(i)));
   
     % calculate orbital period.
    Td = 2*pi*sqrt(ad(i+1)^3/param.mu);

     step = 5*Td; 
    
    tcurr = t; 
    tnext = t+ step;
    
    RAAN(i+1) =  RAAN(i) -param.k*Vdrift^7*cos(Id)*(tnext - tcurr);
    
    t =tnext;
    i = i+1;
    wtg = [wtg,t];
end



if ad(end)-param.Re < 200e3
    'Orbit decayed: choose different drift orbit';
    decaystatus =  true;
else
    decaystatus =  false;
end

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

if isnan(res)
    targetRAAN
    Omega_tf
end

end
