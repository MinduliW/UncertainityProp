function [x, Total_dv1,TOF,Omegadiff ,leg1,leg2, waitTimevec, targetRAAN, wait,fuelburnt]...
    =RAANmatchingmethod(v0,vf, inc0, incf, param)

a0 = param.mu/v0^2;
af = param.mu/vf^2;
x0 = [v0-0.1, inc0-0.001,0.01];

lb = [min([v0,vf])-1000, min([inc0,incf])-20*pi/180, 0 ];
ub = [max([v0,vf]), max([inc0,incf])+20*pi/180, inf];

if param.Topt == false
    ub(end) = 50*param.TU; % kinda need to limit this to be reasonable othewise the guidance will take eons to run
end

[x,~] = fmincon( @(x) fitnessfunction(x, param),x0,[],[],[],[],lb,ub,...
    @(x) condition2(x,param) , optimoptions('fmincon','Display','iter'));


 [~,TOF, Tdv] = fitnessfunction(x, param);

[Total_dv1,TOF,Omegadiff ,leg1,leg2, waitTimevec, targetRAAN, wait,fuelburnt]   = ...
    outcomes(x,0, TOF, param);


% plot result.

figure;
subplot(3,1,1);hold on;
tSeg = linspace(0,TOF, 500);
tT = tSeg*param.TU/86400;
tL1 = leg1.t*param.TU/86400;
wt = leg1.t(end)*param.TU/86400; +wait.timevec*param.TU/86400;
tL2 = leg2.t*param.TU/86400;

plot(tT,  (a0-param.Re)*param.LU/1e3*ones(size(tSeg)), 'g--');
plot(tT, (af-param.Re)*param.LU/1e3*ones(size(tSeg)),'--', 'color', '#C55B11');
plot(tL1,(leg1.a-param.Re)*param.LU/1e3, 'k');
plot(wt,(leg1.a(end)-param.Re)*param.LU/1e3*ones(size(wt)), 'k' )
p = plot(tL2,(leg2.a-param.Re)*param.LU/1e3, 'k');

plot_latex(p, 'time(days)', 'a (km)','','' ,{'Initial','Target','Thrust phase 1',...
    'Coast phase',...
    'Thrust phase 2'});

subplot(3,1,2); hold on;
plot(tT, rad2deg(param.x0(2)*ones(size(tSeg))), 'g--');
plot(tT, rad2deg(param.xf(2)*ones(size(tSeg))),'--', 'color', '#C55B11');
plot(tL1,leg1.inc*180/pi ,'k');
plot(wt,leg1.inc(end)*180/pi*ones(size(wt)) ,'k')
p = plot(tL2,leg2.inc*180/pi, 'k');

plot_latex(p, 'time(days)', 'inc(deg)','', '' ,{ });

Omega_t2 = wrapTo2Pi(leg1.RAAN(end)-param.k*x(1)^7*cos(x(2))*(wait.timevec)); %rad
Omega_x0 = wrapTo2Pi(param.x0(3)-param.k*v0^7*cos(inc0)*tSeg);
Omega_xf =  wrapTo2Pi(param.xf(3) -param.k*param.xf(1)^7*cos(param.xf(2))*(tSeg));

subplot(3,1,3); hold on;
plot(tT, rad2deg(Omega_x0), 'g--');
plot(tT, rad2deg(Omega_xf),'--', 'color', '#C55B11')
plot(tL1,(leg1.RAAN)*180/pi,'k')
plot(wt,Omega_t2*180/pi, 'k');
p = plot(tL2,wrapTo2Pi(leg2.RAAN)*180/pi, 'k');
plot_latex(p, 'time(days)', 'RAAN(deg)','', '',{});

    
    
fprintf('Edelbaum \n')
fprintf('Fuel burnt: %f kg \n', fuelburnt*param.MU);
fprintf('dv: %f m/s \n', Total_dv1*param.LU/param.TU);
fprintf('TOF (d) : %f  \n', TOF*param.TU/86400);
fprintf('--------------------------------- \n')
end
