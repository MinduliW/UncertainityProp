function plottrajwaiting(a0, af, inc0, x, TOF1 ,leg1,leg2,coemean, ...
    times,waitdata,param)

waita = waitdata.a;
v0 = sqrt(param.mu/a0);
vf = sqrt(param.mu/af);

figure;
subplot(4,1,1);hold on;
wt = linspace(leg1.t(end), leg2.t(1), 100);
tSeg = linspace(0,TOF1, 500);

plot(tSeg/param.TU, (a0-param.Re)/1e3*ones(size(tSeg)), 'r--');
plot(tSeg/param.TU, (af-param.Re)/1e3*ones(size(tSeg)),'b--');
plot(times/param.TU, (coemean(:,1)-param.Re)/1e3,'m');
plot(leg1.t/param.TU,(leg1.a-param.Re)/1e3, 'g--');

if param.drag == false
    plot(wt/param.TU,(leg1.a(end)-param.Re)/1e3*ones(size(wt)), 'g--' )
else
    
    plot(wt/param.TU,(waita-param.Re)/1e3, 'g--' )
end

    
p = plot(leg2.t/param.TU,(leg2.a-param.Re)/1e3, 'g--');



plot_latex(p, 'time(days)', 'a (km)','', '' ,{'Initial','Target','Forward Prop.',...
    'Edelbaum'});

subplot(4,1,2); hold on;
plot(tSeg/param.TU, rad2deg(param.x0(2)*ones(size(tSeg))), 'r--');
p = plot(tSeg/param.TU, rad2deg(param.xf(2)*ones(size(tSeg))),'b--');
p = plot(times/param.TU, coemean(:,3)*180/pi,'m');
plot(leg1.t/param.TU,leg1.inc*180/pi ,'g--');
plot(wt/param.TU,leg1.inc(end)*180/pi*ones(size(wt)) ,'g--')
p = plot(leg2.t/param.TU,leg2.inc*180/pi, 'g--');


plot_latex(p, 'time(days)', 'inc(deg)','', '' ,{ });

Omega_t2 = unwrap(leg1.RAAN(end) -param.k*x(1)^7*cos(x(2))*(wt-leg1.t(end))); %rad
Omega_x0 = unwrap(param.x0(3)-param.k*v0^7*cos(inc0)*tSeg);
Omega_xf =  unwrap(param.xf(3) -param.k*param.xf(1)^7*cos(param.xf(2))*(tSeg));

subplot(4,1,3); hold on;
plot(tSeg/param.TU, rad2deg(Omega_x0), 'r--');
plot(tSeg/param.TU, rad2deg(Omega_xf),'b--')
p = plot(times/param.TU, unwrap(coemean(:,4))*180/pi,'m');
plot(leg1.t/param.TU,unwrap(leg1.RAAN)*180/pi,'g--')
plot(wt/param.TU,Omega_t2*180/pi, 'g--');
p = plot(leg2.t/param.TU,unwrap(leg2.RAAN)*180/pi, 'g--');

plot_latex(p, 'time(days)', '$\Omega$ (deg)','', '' ,...
    {});


subplot(4,1,4); hold on; 

plot(times/param.TU, coemean(:,2),'m');
plot_latex(p ,'time (days)', 'eccentricity', '', '', {})




end
