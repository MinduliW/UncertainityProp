function [result,totalSmaChange, totalIncChange,Total_TOF, Total_dV,...
    Total_fuel,TotalThrustTime,decaystatus] = tourGeneral(x, debris,m_debris, param)

totalSmaChange = 0;
totalIncChange = 0;
%param.t0 = x(end); 
datevec = '19-June-2024';
param.t0 = juliandate(datevec,'dd-mmm-yyyy'); %+ 1000;

param.timeprev = 0; 
param.legno = 1; 
% Leg 1 : From debris orbit 1 to 350 km orbit, No RAAN matching required. 
param.m0  = param.m_wet_servicer + m_debris(1);

% from debris orbit 
a0 = debris(1,1);inc0 = debris(1,2);RAAN0 = debris(1,3); 

% to deorbit orbit
af = param.target_altitude + param.Re; incf = inc0;

totalSmaChange = totalSmaChange + abs(af - a0);
totalIncChange = totalIncChange + abs(incf -inc0);


[TOF,Tdv, lastOmega,semiMajor, inclination, RAAN,tSeg,~,~,~,mass]=kluver(a0,  af, inc0,incf,RAAN0,param);

TOFs = [TOF/param.TU];
dvs = [Tdv(end)/1e3]; 

Total_TOF = TOF/param.TU; Total_dV = Tdv(end)/1e3; 


TotalThrustTime = TOF/param.TU*param.dutyRatio; 


mass_after = param.m0/exp(Tdv(end)/param.Isp/param.g0); 

fuelexp = [param.m0-mass_after];
wt = zeros(1,length(m_debris)-1);

if param.plots == true
    figure;
    subplot(4,1,1);hold on;
    plot(tSeg/param.TU, a0/1e3*ones(size(tSeg)), '--');
    plot(tSeg/param.TU, af/1e3*ones(size(tSeg)),'--');
    p = plot(tSeg/param.TU,semiMajor/1e3);
    plot_latex(p, 'time(days)', 'a (km)','', 'Leg1 : From Debris 1 to orbit below ISS' ,{'$a_0$','$a_f$','$a(t)$'});
    
    subplot(4,1,2); hold on;
    plot(tSeg/param.TU, rad2deg(inc0*ones(size(tSeg))), '--');
    plot(tSeg/param.TU, rad2deg(incf*ones(size(tSeg))),'--');
    p = plot(tSeg/param.TU,inclination*180/pi);
    plot_latex(p, 'time(days)', 'inc(deg)','', '' ,{'$i_0$','$i_f$','$i(t)$'});
    
    subplot(4,1,3); hold on;
    Omega_x0 = RAAN0 -param.k*(param.mu/a0)^3.5*cos(inc0)*tSeg;
    plot(tSeg/param.TU, rad2deg(Omega_x0), '--');
    p = plot(tSeg/param.TU,RAAN*180/pi);
    plot_latex(p, 'time(days)', '$\Omega$(deg)','', '' ,{'$\Omega_0$','$\Omega(t)$'});
    
    subplot(4,1,4); hold on;
    p = plot(tSeg/param.TU,mass);
    plot_latex(p, 'time(days)', 'mass (kg)','', '' ,{});
    
    
end



k = 1; 
currstart = 1; 
for j = 1: length(m_debris)-1
    %% Proximity operations to handover the debris to the sheperd.
    Omega_0_new = lastOmega -param.k*(param.mu/af)^3.5*cos(incf)*param.RDV2;
    
    Total_TOF = Total_TOF + param.RDV2/param.TU;
   
    TOFs = [TOFs,param.RDV2/param.TU];
    dvs = [dvs, 0];

    %% From deorbit orbit to debris orbit, RAAN matching required.
     param.timeprev = Total_TOF; 
    param.legno = param.legno + 1; 
    param.m0 = mass_after -  m_debris(j);
    
    % from deorbit orbit
    a0 = param.target_altitude + param.Re; inc0 = incf; RAAN0 = Omega_0_new;
    
    % to debris orbit 2
    af =  debris(j+1,1);
    incf = debris(j+1,2);
%    RAANf_t0 = Omega_0_new + param.dRAAN(j);
    RAANf_t0 = debris(j+1,3)-param.k*(param.mu/af)^3.5*cos(incf)*Total_TOF*param.TU;
    if param.optimiseRAAN == true
      RAANf_t0 = (x(3*j))-param.k*(param.mu/af)^3.5*cos(incf)*Total_TOF*param.TU;
    end
    
        
    totalSmaChange = totalSmaChange + abs(af - a0);
    totalIncChange = totalIncChange + abs(incf -inc0);
    
    if param.optimiseRAAN == true
        
        [Tdv,TOF,Omega2,wt(k),~,thrusttime,decaystatus] = propOpt(x(currstart:currstart+1), a0,inc0,RAAN0 , af, incf,...
            RAANf_t0, param);
        currstart = currstart + 3;
    else
        [Tdv,TOF,Omega2,wt(k),~,thrusttime,decaystatus] = propOpt(x(j+k -1:j+k), a0,inc0,RAAN0 , af, incf,...
            RAANf_t0, param);
    end
    
    TotalThrustTime = TotalThrustTime+thrusttime;
    k = k +1; 
    
    Total_TOF = Total_TOF +TOF/param.TU;
    Total_dV = Total_dV + Tdv(end)/1e3; 
    
    TOFs = [TOFs,TOF/param.TU];
    dvs = [dvs, Tdv(end)/1e3];
      
    mass_after = param.m0 /exp(Tdv(end)/param.Isp/param.g0);

    fuelexp = [fuelexp, param.m0-mass_after];

    %% Proximity operations at the target
    RAAN2 = Omega2 -param.k*(param.mu/af)^3.5*cos(incf)*param.RDV1;
    Total_TOF = Total_TOF +param.RDV1/param.TU;
   
     TOFs = [TOFs,param.RDV1/param.TU];
     dvs = [dvs, 0];
      
    
    % Leg 3 
    param.timeprev = Total_TOF; 
    param.legno = param.legno + 1; 
    param.m0 = mass_after +  m_debris(j+1);
    % from debris orbit 2
    a0 = debris(j+1,1); inc0 = debris(j+1,2); RAAN0 = RAAN2;
    
    % to deorbit orbit
    af = param.target_altitude + param.Re; incf = inc0;
    
    totalSmaChange = totalSmaChange + abs(af - a0);
    totalIncChange = totalIncChange + abs(incf -inc0);

    
    [TOF,Tdv,lastOmega,semiMajor, inclination, RAAN,tSeg,~,~,~,mass] =  kluver(a0,  af, inc0,incf,RAAN0,param);
    
    Total_TOF = Total_TOF +TOF/param.TU;
    Total_dV = Total_dV + Tdv(end)/1e3; 
  
    TOFs = [TOFs,TOF/param.TU];
    dvs = [dvs, Tdv(end)/1e3];
    
    
    if param.plots==true
        figure;
        subplot(4,1,1);hold on;
        plot(tSeg/param.TU, a0/1e3*ones(size(tSeg)), '--');
        plot(tSeg/param.TU, af/1e3*ones(size(tSeg)),'--');
        p = plot(tSeg/param.TU,semiMajor/1e3);
        
        plot_latex(p, 'time(days)', 'a (km)','',  strcat('Leg' , '{ }', ...
            num2str(param.legno), ': From Debris to orbit below ISS') ,{'$a_0$','$a_f$','$a(t)$'});
        
        subplot(4,1,2); hold on;
        plot(tSeg/param.TU, rad2deg(inc0*ones(size(tSeg))), '--');
        plot(tSeg/param.TU, rad2deg(incf*ones(size(tSeg))),'--');
        p = plot(tSeg/param.TU,inclination*180/pi);
        plot_latex(p, 'time(days)', 'inc (deg)','', '' ,{'$i_0$','$i_f$','$i(t)$'});
        
        subplot(4,1,3); hold on;
        Omega_x0 = RAAN0 -param.k*(param.mu/a0)^3.5*cos(inc0)*tSeg;
        plot(tSeg/param.TU, rad2deg(Omega_x0), '--');
        p = plot(tSeg/param.TU,RAAN*180/pi);
        plot_latex(p, 'time(days)', '$\Omega$(deg)','', '' ,{'$\Omega_0$','$\Omega(t)$'});
        
        subplot(4,1,4); hold on;
        p = plot(tSeg/param.TU,mass);
        plot_latex(p, 'time(days)', 'mass (kg)','', '' ,{});
        
    end
    
    
    TotalThrustTime = TotalThrustTime+TOF/param.TU*param.dutyRatio;
    
    mass_after = param.m0/exp(Tdv(end)/param.Isp/param.g0);
    
    
    fuelexp = [fuelexp, param.m0-mass_after];

end



% Proximity operations to handover the debris to the sheperd.
Total_TOF = Total_TOF +param.RDV2/param.TU;

% TOFs = [TOFs,param.RDV2/param.TU];
% dvs = [dvs, 0];
Total_fuel = sum(fuelexp);
% TOFs = [TOFs';Total_TOF]
% 
% 
% dvs = [dvs'*1e3;Total_dV*1e3]
% 
%   
if param.Topt == true
%   result = 2*Total_TOF*stepfcn(Total_dV, param.dvLimit)*(Total_dV - param.dvLimit)^2 ...
%       + Total_TOF;
  result = Total_TOF; 
else
%    result = stepfcn(max(wt), param.waitTimeLimit)*(max(wt) -...
%         param.waitTimeLimit)^2 + Total_dV;  
   result = Total_fuel; %Total_dV;
end

if decaystatus > 0
    result = 1e10;
end



