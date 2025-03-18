function [TOF,wt] = getTOF(x, k, param)

Vd = x(1);
Id = x(2);

a0 = param.mu/param.x0(1)^2; %m
ad = param.mu/Vd^2; %m
af = param.mu/param.xf(1)^2; %m

% Calculate the RAAN change during the first burn
if Vd ~= param.x0(1)
    [t1,~, deltaOmega1] =  kluver(a0, ad, param.x0(2),Id,0,param);
else
    t1 = 0; deltaOmega1 = 0;
end

% RAAN derivative of waiting 
Omega_wait_dot =  -param.k*Vd^7*cos(Id);

% Calcualte the RAAN change during second burn 
if Vd ~= param.xf(1)
   [t2,~, deltaOmega2] =  kluver(ad,af, Id,param.xf(2) ,0,param);
else
    t2 = 0; deltaOmega2 = 0;
end

% RAAN derivative of the target 
Omega_target_dot = -param.k*param.xf(1)^7*cos(param.xf(2)); 

wt = (param.x0(3) + deltaOmega1 +deltaOmega2 - param.xf(3) - ...
    Omega_target_dot*(t1 + t2)- 2*k*pi)/(Omega_target_dot - Omega_wait_dot);


if wt <0
    wt = 100000*param.TU;
end

TOF = wt +t1 + t2; 


end
