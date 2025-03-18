function [result,TOF, Tdv] = fitnessfunction(x, param)

Vd = x(1);
Id = x(2);
waitTime = x(3); 

a0 = param.mu/param.x0(1)^2; %m
ad = param.mu/Vd^2; %m
af = param.mu/param.xf(1)^2; %m


%% Leg 1 
param.timeprev = 0; 
[t1,Tdv1, Omega_t1,~, ~, ~,~,~,~,~,massafterl1] =  kluver(a0, ad, param.x0(2),Id,...
    param.x0(3),param);

%% Drifting
[Omega_t2, semiMajor, dv_dragCorrection, masscons] = waitDragEffect(Omega_t1, ad, Id,waitTime,massafterl1(end), param);

%% Leg 2 
param.m0 = massafterl1(end)-masscons;
param.timeprev = t1 + waitTime; 
[t2,Tdv2, ~] =  kluver(semiMajor(end), af, Id,param.xf(2) ,...
   Omega_t2,param);

%% Totals 
TOF = t1+ t2 + waitTime;

Tdv = Tdv1(end) + Tdv2(end)+dv_dragCorrection; 

if param.Topt == true
result = TOF; 
else
 result = Tdv; 
end

end


