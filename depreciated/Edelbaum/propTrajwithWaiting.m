function [times, states,dv,fuelburnt,tof] = propTrajwithWaiting(a0, RAAN0_t0, inc0, omega0, e0,... ...
    incf, leg1, leg2,waitTimevec,x, param)

isguidance = param.guidance; 

%% Leg 1 
mintheta = 0;
diffold = 100;
for t = 0:0.0001:2*pi
    diff = getx0(t, a0, RAAN0_t0, inc0,param);
    
    if abs(diff) < abs(diffold)
        mintheta = t;
    end
    
    diffold = diff;
end

x0 = CoordConv.kepler2MEOE([a0, e0,omega0,RAAN0_t0,inc0,mintheta]);
period =  2*pi*sqrt(a0^3/param.mu);

param.Thrust = true;

if x(2) > inc0
    param.incincrease = true;
else
    param.incincrease = false;
end

param.Edelsma = leg1.a;
param.Edelinc = leg1.inc;
param.EdelRAAN = leg1.RAAN;

if isguidance == true
    param.guidance = true;
end

opt   = odeset('Events', @(t,x) myEvent(t,x,param.Edelsma(end),param));
if abs(leg1.t(1)- leg1.t(end)) >100
    [time_fp_leg1,states_fp_leg1] = ode45(@(t,x) forwardPropagator(t,x, ...
        leg1.beta,leg1.fs, leg1.t, param),leg1.t(1): period/20: leg1.t(end),...
        [x0';param.m0],opt);
else
    time_fp_leg1 = leg1.t(1);
    states_fp_leg1 = [x0';param.m0]';
end

fprintf('Errors after Leg 1 \n')
errorcalc(time_fp_leg1, states_fp_leg1, param.Edelsma(end),param.Edelinc(end),...
    param.EdelRAAN(end),param);
fprintf('___________________________\n')

dv1 = param.Isp*param.g0*log(param.m0/states_fp_leg1(end,end));



%% Waiting 

param.Edelsma = leg1.a(end)*ones(size(waitTimevec));
param.Edelinc = leg1.inc(end)*ones(size(waitTimevec));
param.EdelRAAN =leg1.RAAN(end) -param.k*x(1)^7*cos(x(2))*waitTimevec; %rad


if abs(waitTimevec(1)- waitTimevec(end)) >100
    param.Thrust = false;
%     if param.guidance == true
%         param.guidance = false;
%     end
    
    
    opt   = odeset('Events', @(t,x) myEvent2(t,x,param.EdelRAAN(end),param));
    
    [time_fp_wt,states_fp_wt] = ode45(@(t,x) forwardPropagator(t,x, [],[],...
        waitTimevec, param),[waitTimevec(1), waitTimevec(end)],...
        states_fp_leg1(end,:),opt);
else
    time_fp_wt = waitTimevec(1);
    states_fp_wt = states_fp_leg1(end,:);
end


fprintf('Errors after waiting \n')
errorcalc(time_fp_wt, states_fp_wt, param.Edelsma(end),param.Edelinc(end),...
    param.EdelRAAN(end),param);
fprintf('___________________________\n')


% times = [time_fp_leg1;time_fp_wt];
% states = [states_fp_leg1;states_fp_wt];
% return;


%% Leg 2

if  incf> x(2)
    param.incincrease = true;
else
    param.incincrease = false;
end


param.Edelsma = leg2.a;
param.Edelinc = leg2.inc;
param.EdelRAAN =leg2.RAAN;


param.Thrust = true;
if isguidance == true
    param.guidance = true;
end


opt   = odeset('Events', @(t,x) myEvent(t,x,param.Edelsma(end),param));
[time_fp_leg2,states_fp_leg2] = ode45(@(t,x) forwardPropagator...
    (t,x, leg2.beta,leg2.fs, leg2.t, param),...
    [leg2.t(1): period/20:leg2.t(end)],[states_fp_wt(end,:)],opt);


dv2 = param.Isp*param.g0*log(states_fp_wt(end,end)/states_fp_leg2(end,end));

fprintf('Errors after Leg2 \n')
errorcalc(time_fp_leg2, states_fp_leg2, param.Edelsma(end),param.Edelinc(end),...
    param.EdelRAAN(end),param);
fprintf('___________________________\n')


%% Totals 

times = [time_fp_leg1;time_fp_wt; time_fp_leg2];
states = [states_fp_leg1;states_fp_wt; states_fp_leg2];

dv = dv1 + dv2;
tof = (time_fp_leg2(end) - time_fp_leg1(1))/param.TU;

fuelburnt = param.m0 - states_fp_leg2(end,end);

end

