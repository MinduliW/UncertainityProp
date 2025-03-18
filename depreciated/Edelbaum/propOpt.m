function [Total_dv1,TOF1,targetRAAN,wT,fuelconsumed,thrusttime] = propOpt...
    (x , a0,inc0,RAAN0 , af, incf, RAANf,param)


% Initital orbit at the start of leg 2
v0  =  sqrt(param.mu/a0);
param.x0 = [v0 , inc0,RAAN0];

% Target orbit at the start of leg 2
vf  = sqrt(param.mu/af);
param.xf = [vf , incf, RAANf];


ks = param.krange;
wts = zeros(size(ks));
for i = 1: length(ks)
    % calculate waittime in days
    [wts(i)] = initWaitTime2(x,ks(i), param);
end

%decaystatus
% if all(decaystatus>0)
%     Total_dv1 = 0; TOF1 = 0;targetRAAN =0;waittime = 0;rrs =0;thrusttime=0;
%     decaystatus =1;
%     return; 
% end



waittimes = wts(wts>0);
ks = ks(wts>0);

[waittime, index] = min(waittimes);
k = ks(index);



x(3) =  waittime;
% if param.addDrag == true
%     a0 = param.mu/param.x0(1)^2; %m
%     ad = param.mu/x(1)^2; %m
%     af = param.mu/param.xf(1)^2; %m
%     
%     % Calculate the RAAN change during the first burn
%     [t1,~, Omega1,~,~,~,~,~,~,~,massL1] =  kluver(a0, ad, param.x0(2),x(2),...
%         param.x0(3),param);
%     
%     
%     %waittime = sqwaittimenodrag^2*param.TU;
%     
%     % calculate the thrust required to maintain the altitude
%     
%     
%     options = optimoptions('fsolve', 'display', 'none');
%     [sqwaittime,~,~] = fsolve(@(wt) solveWithDrag(wt, x(1),x(2), Omega1,...
%         massL1(end),t1(end), param),...
%         sqrt(waittime/param.TU),options);
%     
%     x(3) = (sqwaittime^2*param.TU);
% else
%     x(3) =  waittime;
% end

wT = x(3);

[Total_dv1,TOF1,~ ,leg1,leg2,~, targetRAAN,wait] = ...
    outcomes(x,k , param);



thrusttime = leg1.t1/param.TU*param.dutyRatio +  leg2.t2/param.TU*param.dutyRatio;
fuelconsumed = leg1.mass(1)-leg2.mass(end);


if param.plots == true
    
     param.timeprev  =  param.TOFbefore;
   
    figure;
    subplot(3,1,1);hold on;
    %wt = linspace(leg1.t(end), leg2.t(1), 500)/param.TU;
    tSeg = linspace(0,TOF1, 500);
    tT = param.timeprev + tSeg/param.TU; 
    tL1 = param.timeprev +leg1.t/param.TU;
    wt = param.timeprev+ leg1.t(end)/param.TU +wait.timevec/param.TU;
    tL2 = param.timeprev +leg2.t/param.TU;
    
    param.timeprev  = round(param.timeprev );
    plot(tT,  (a0-param.Re)/1e3*ones(size(tSeg)), 'g--');
    plot(tT, (af-param.Re)/1e3*ones(size(tSeg)),'--', 'color', '#C55B11');
    plot(tL1,(leg1.a-param.Re)/1e3, 'k');
    plot(wt,(leg1.a(end)-param.Re)/1e3*ones(size(wt)), 'k' )
    p = plot(tL2,(leg2.a-param.Re)/1e3, 'k');
    
    if param.guidance == true
        p = plot(timeG/param.TU+param.timeprev ,coemean(:,1)/1e3, 'm-');
    end
    xlim([min(param.timeprev +tSeg/param.TU),max(param.timeprev +tSeg/param.TU)]);
   
    plot_latex(p, 'time(days)', 'a (km)','', strcat('Leg' , '{ }', ...
        num2str(param.legno), ': From Shepherd orbit to debris orbit') ,{'Initial','Target','Thrust phase 1',...
        'Coast phase',...
        'Thrust phase 2', 'Guidance'});
    
    subplot(3,1,2); hold on;
    plot(tT, rad2deg(param.x0(2)*ones(size(tSeg))), 'g--');
    plot(tT, rad2deg(param.xf(2)*ones(size(tSeg))),'--', 'color', '#C55B11');
    plot(tL1,leg1.inc*180/pi ,'k');
    plot(wt,leg1.inc(end)*180/pi*ones(size(wt)) ,'k')
    p = plot(tL2,leg2.inc*180/pi, 'k');
    if param.guidance == true
        p = plot(param.timeprev +timeG/param.TU,coemean(:,3)*180/pi, 'm-');
    end
    %ylim([98,99])
    xlim([min(tT),max(tT)]);
   
    plot_latex(p, 'time(days)', 'inc(deg)','', '' ,{ });
    

 
    Omega_t2 = wrapTo2Pi(leg1.RAAN(end)-param.k*x(1)^7*cos(x(2))*(wait.timevec)); %rad
    Omega_x0 = wrapTo2Pi(param.x0(3)-param.k*v0^7*cos(inc0)*tSeg);
    Omega_xf =  wrapTo2Pi(param.xf(3) -param.k*param.xf(1)^7*cos(param.xf(2))*(tSeg) + 2*k*pi);
    
    subplot(3,1,3); hold on;
    plot(tT, rad2deg(Omega_x0), 'g--');
    plot(tT, rad2deg(Omega_xf),'--', 'color', '#C55B11')
    plot(tL1,(leg1.RAAN)*180/pi,'k')
    plot(wt,Omega_t2*180/pi, 'k');
    p = plot(tL2,wrapTo2Pi(leg2.RAAN)*180/pi, 'k');
    if param.guidance == true
        p = plot(param.timeprev +timeG/param.TU,wrapTo2Pi(coemean(:,4))*180/pi, 'm-');
    end
    ylim([0,360]);
    xlim([min(param.timeprev +tSeg/param.TU),max(param.timeprev +tSeg/param.TU)]);
   
    plot_latex(p, 'time(days)', '\Omega(deg)','', strcat('Leg' , '{ }', ...
        num2str(param.legno), ': From orbit below ISS to debris orbit') ,{'Initial','Target','Thrust phase 1',...
        'Coast phase',...
        'Thrust phase 2', 'Guidance'});

    
%     subplot(4,1,4); hold on;
%     plot(tL1,leg1.mass,'k');
%     plot(wt,leg1.mass(end)*ones(size(wt)),'k');
%     plot(tL2,leg2.mass, 'k');
%     if param.guidance == true
%     p = plot(param.timeprev +timeG/param.TU,coemean(:,7), 'm-');
%     end
%   
%     xlim([min(param.timeprev +tSeg/param.TU),max(param.timeprev +tSeg/param.TU)]);
%    
%     plot_latex(p, 'time(days)', 'Mass_{serv.}(kg)','', '' ,{});
%     
%     subplot(5,1,5); hold on;
%     plot(param.timeprev +leg1.t/param.TU,param.dVbefore*1e3 + leg1.dV,'k');
%     plot(param.timeprev +wt/param.TU,param.dVbefore*1e3 + leg1.dV(end)*ones(size(wt)),'k');
%     p = plot(param.timeprev +leg2.t/param.TU,param.dVbefore*1e3 + leg1.dV(end) + leg2.dV,'k');
%     if param.guidance == true
%         DVg = param.Isp*param.g0*log(param.m0./coemean(:,7));
%         plot(param.timeprev +timeG/param.TU,param.dVbefore*1e3 + DVg, 'm-');
%     end
%     
%     xlim([min(param.timeprev +tSeg/param.TU),max(param.timeprev +tSeg/param.TU)]);
%    
   % plot_latex(p, 'time(days)', '\Delta v (m/s)','', '' ,{});
    
    


end

end