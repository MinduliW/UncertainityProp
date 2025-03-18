function [c, ceq]  = condition2(x, param)
Vd = x(1);
Id = x(2);
waitTime = x(3); 

a0 = param.mu/param.x0(1)^2; %m
ad = param.mu/Vd^2; %m
af = param.mu/param.xf(1)^2; %m

%% Leg 1 
param.timeprev = 0; 
[t1,~, Omega_t1,~, ~, ~,~,~,~,~,massafterl1] =  kluver(a0, ad, param.x0(2),Id,...
    param.x0(3),param);

%% Waiting 
Omega_t2 = Omega_t1 -param.k*Vd^7*cos(Id)*(waitTime);

%% Leg 2 
param.m0 = massafterl1(end);
param.timeprev = t1 + waitTime; 
[t2,~, Omega_tf] =  kluver(ad, af, Id,param.xf(2) ,...
   Omega_t2,param);

%% Omega evolution.
Omega_target_dot = -param.k*param.xf(1)^7*cos(param.xf(2));

RAANf_tf = param.xf(3)+ Omega_target_dot*( t1+ t2 + waitTime);

ceq = sin(0.5*(Omega_tf - RAANf_tf));
c = [];
end

