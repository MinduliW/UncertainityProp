
function ceq = solveforRAAN(x, Vd, Id, param)

waitTime = x*param.TU;

a0 = param.mu/param.x0(1)^2; %m
ad = param.mu/Vd^2; %m
af = param.mu/param.xf(1)^2; %m


% Calculate the RAAN change during the first burn
if Vd ~= param.x0(1)
    [t1,~, Omega_t1] =  kluver(a0, ad, param.x0(2),Id,...
        param.x0(3),param);
else
    t1 = 0; Omega_t1 =  param.x0(3);
end


Omega_t2 = Omega_t1 -param.k*Vd^7*cos(Id)*(waitTime);


if Vd ~= param.xf(1)
    
    [t2,~, Omega_tf] =  kluver(ad, af, Id,param.xf(2) ,...
        Omega_t2,param);
    
else
    t2 = 0; Omega_tf = Omega_t2;
end


finalTime = t1+ t2 + waitTime;

% evolve final omega.
Omega_target_dot = -param.k*param.xf(1)^7*cos(param.xf(2));

RAANf_tf = param.xf(3)+ Omega_target_dot*finalTime; % correct


difference = (Omega_tf) - (RAANf_tf);
ceq = sin(0.5*(difference)); % difference*180/pi; %;

end