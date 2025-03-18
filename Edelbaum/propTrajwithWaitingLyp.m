function [times, states,dv,fuelburnt,tRaan,TOF] = propTrajwithWaitingLyp(a0, RAAN0_t0, inc0, ...
    incf,  omega0, e0, leg1, leg2,waitTimevec,x, param)
addpath('/Users/minduli/Astroscale_ADR/Main/Applications/Guidance/SlimLocoche');


%% Leg 1 
mintheta = 0;
diffold = 100;
 param.tracktarget =false;
param.followEdelbaum = true;
param.ef = 0;

for t = 0:0.0001:2*pi
    diff = getx0(t, a0, RAAN0_t0, inc0,param);
    
    if abs(diff) < abs(diffold)
        mintheta = t;
    end
    
    diffold = diff;
end

x0 = CoordConv.kepler2MEOE([a0, e0,omega0,RAAN0_t0,inc0,mintheta]);

param.Thrust = true;

param.Edelsma = leg1.a;
param.Edelinc = leg1.inc;
param.EdelRAAN = leg1.RAAN;

param.af = param.Edelsma(end); 
param.incf = param.Edelinc(end); 
param.RAANf = param.EdelRAAN(end);

param.ef = 0; 


if abs(leg1.t(1)- leg1.t(end)) >100
    
    cf = [0.01, 0.01, 1, 0, 0, 1];
    opt   = odeset('Events', @(t,y) stopCriterion(t,y,cf,leg1.t, param));
    
    [time_fp_leg1,states_fp_leg1] = ode45(@(t,x) propagate(t,x,cf(1), cf(2),...
        cf(3), cf(4), cf(5), cf(6),leg1.t,param),[0, leg1.t(end)],[x0';param.m0],opt);
else
    time_fp_leg1 = leg1.t(1);
    states_fp_leg1 = [x0';param.m0]';
end


fprintf('Errors after Leg 1 \n')

tRaan = param.EdelRAAN(end); 
coemean = errorcalc(time_fp_leg1, states_fp_leg1, param.Edelsma(end),param.Edelinc(end),...
    tRaan,param);
fprintf('___________________________\n')

dv1 = param.Isp*param.g0*log(param.m0/states_fp_leg1(end,end));


%% Waiting 

param.followEdelbaum = false;

param.Edelsma = leg1.a(end)*ones(size(waitTimevec));
param.Edelinc = leg1.inc(end)*ones(size(waitTimevec));
param.EdelRAAN =tRaan -param.k*x(1)^7*cos(x(2))*(waitTimevec-leg1.t(end)); %rad

param.af = param.Edelsma(end); 
param.incf = param.Edelinc(end); 
param.RAANf = param.EdelRAAN(end);
param.ef = 0; 

if abs(waitTimevec(1)- waitTimevec(end)) >100
    param.Thrust = false;
    
    
    %opt   = odeset('Events', @(t,y) stopCriterion(t,y,cf,waitTimevec, param));
    
    [time_fp_wt,states_fp_wt] = ode45(@(t,x) propagate(t,x,cf(1), cf(2),...
        cf(3), cf(4), cf(5), cf(6),waitTimevec,param),[waitTimevec(1),...
        waitTimevec(end)],[states_fp_leg1(end,:)]);

else
    time_fp_wt = waitTimevec(1);
    states_fp_wt = states_fp_leg1(end,:);
end


fprintf('Errors after waiting \n')
tRaan = param.EdelRAAN(end);

coemean = errorcalc(time_fp_wt, states_fp_wt, param.Edelsma(end),param.Edelinc(end),...
    tRaan,param);
fprintf('___________________________\n')


% times = [time_fp_leg1;time_fp_wt];
% states = [states_fp_leg1;states_fp_wt];
% return;


%% Leg 2

param.followEdelbaum = true;
param.Thrust = true;

% cf = [0.01, 0.01, 1, 0, 0, 1];
param.Edelsma = leg2.a;
param.Edelinc = leg2.inc;
param.EdelRAAN =leg2.RAAN;

param.af = param.Edelsma(end); 
param.incf = param.Edelinc(end);

tof = waitTimevec(end) - time_fp_leg1(1); 
param.RAANf = param.EdelRAAN(end); %unwrap(param.xf(3) -param.k*param.xf(1)^7*cos(param.xf(2))*(tof));

opt   = odeset('Events', @(t,y) stopCriterion(t,y,cf,leg2.t, param));
[time_fp_leg2,states_fp_leg2] = ode45(@(t,x) propagate(t,x,cf(1), cf(2),...
    cf(3), cf(4), cf(5), cf(6),leg2.t,param),[leg2.t(1), leg2.t(end)]...
    ,[states_fp_wt(end,:)]);

%+100*24*60*60
stopCriterion(time_fp_leg2(end),states_fp_leg2(end,:),cf,leg2.t, param)

for t = 1: length(time_fp_leg2)
    [~,~,~, L(t)] = stopCriterion(time_fp_leg2(t),...
        states_fp_leg2(t,:),cf,leg2.t, param);
end

figure;
plot(L);

dv2 = param.Isp*param.g0*log(states_fp_wt(end,end)/states_fp_leg2(end,end));

fprintf('Errors after Leg2 \n');
n = sqrt(param.mu./param.af.^3);  %1/s
RAANdot = -3/2*param.J2*n.*(param.Re./param.af).^2.*cos(param.incf);
tof = (time_fp_leg2(end) - time_fp_leg1(1))/param.TU;


%tRaan =param.xf(3) + RAANdot.*(tof);
tRaan=  unwrap(param.xf(3) -param.k*param.xf(1)^7*cos(param.xf(2))*(tof)*param.TU);



errorcalc(time_fp_leg2, states_fp_leg2, param.Edelsma(end),param.Edelinc(end),...
    tRaan,param);
fprintf('___________________________\n')


%% Totals 

times = [time_fp_leg1;time_fp_wt; time_fp_leg2];
states = [states_fp_leg1;states_fp_wt; states_fp_leg2];

dv = dv1 + dv2;

fuelburnt = param.m0 - states_fp_leg2(end,end);
TOF = times(end)/param.TU;
end
