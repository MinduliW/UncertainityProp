function tourGeneralOutcome(x, debris,m_debris, param)

 
datevec = '19-June-2024';
param.t0 = juliandate(datevec,'dd-mmm-yyyy'); %+ 1000;

fprintf('Launch date : %s \n',  datestr(datevec));
fprintf('Isp :  %.2f s  \n', param.Isp); 
fprintf('Maximum thrust : %.2f N \n', param.T); 

if param.eclipses == true
    fprintf('Trajectory accounts for eclipses \n'); 
else
    fprintf('Duty Ratio : %.2f  \n', param.dutyRatio); 
end

fprintf('Target altitude = %.2f km \n', param.target_altitude/1e3);
fprintf('..............................\n')


%% Leg 1 : From debris orbit 1 to 350 km orbit, No RAAN matching required. 
param.m0  = param.m_wet_servicer + m_debris(1);

param.legno = 1; 

% from debris orbit 
a0 = debris(1,1);inc0 = debris(1,2);RAAN0 = debris(1,3); 

% to deorbit orbit
af = param.target_altitude + param.Re; incf = inc0;


[TOF,Tdv, lastOmega,semiMajor, inclination, RAAN,tSeg,~,~,~,mass]=...
    kluver(a0,  af, inc0,incf,RAAN0,param);


Total_TOF = TOF/param.TU; Total_dV = Tdv(end)/1e3; 

mass_after = param.m0/exp(Tdv(end)/param.Isp/param.g0); 

mf(1) = param.m0  - mass_after;


fprintf('Leg 1: From Debris 1 to orbit below ISS \n')
fprintf('From a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', a0/1e3, inc0*180/pi,...
    RAAN0*180/pi);
fprintf('To a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
    lastOmega*180/pi);
fprintf('Delta v = %f km/s \n', Tdv(end)/1e3);
fprintf('TOF = %f days \n', TOF/param.TU);

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
yyaxis right; 
p = plot(tSeg/param.TU,Tdv);
plot_latex(p, 'time(days)', '$\Delta v$ (m/s)','', '' ,{});




wt = zeros(1,length(m_debris)-1);

k = 1; 
currstart = 1; 
for j = 1: length(m_debris)-1
    %% Proximity operations to handover the debris to the sheperd.
    Omega_0_new = lastOmega -param.k*(param.mu/af)^3.5*cos(incf)*param.RDV2;
    
    Total_TOF = Total_TOF + param.RDV2/param.TU;

    fprintf('..............................\n')
    fprintf('Proximity operations to handover the debris %d to the sheperd. \n', i-1)
    fprintf('From a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
        lastOmega*180/pi);
    fprintf('To a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
        Omega_0_new*180/pi);
    fprintf('TOF = %f days \n', param.RDV2/param.TU);


    %% From deorbit orbit to debris orbit, RAAN matching required.
    
    param.legno =  param.legno + 1; 
    param.m0 = mass_after -  m_debris(j);
    
    % from deorbit orbit
    a0 = param.target_altitude + param.Re; inc0 = incf; RAAN0 = Omega_0_new;
    
    % to debris orbit 2
    af =  debris(j+1,1);
    incf = debris(j+1,2);
    RAANf_t0 = debris(j+1,3)-param.k*(param.mu/af)^3.5*cos(incf)*Total_TOF*param.TU;
    if param.optimiseRAAN == true
      RAANf_t0 = x(3*j)-param.k*(param.mu/af)^3.5*cos(incf)*Total_TOF*param.TU;
    end
   
    if param.optimiseRAAN == true
        [Tdv,TOF,Omega2,wt(k)] = propOpt(x(currstart:currstart+1), a0,inc0,RAAN0 , af, incf,...
            RAANf_t0, param);
        currstart = currstart + 3;
    else
        [Tdv,TOF,Omega2,wt(k)] = propOpt(x(j+k -1:j+k), a0,inc0,RAAN0 , af, incf,...
            RAANf_t0, param);
    end
    
    k = k +1; 
    
    Total_TOF = Total_TOF +TOF/param.TU;
    Total_dV = Total_dV + Tdv/1e3; 
    mass_after = param.m0 /exp(Tdv/param.Isp/param.g0);

    mf= [mf, param.m0  - mass_after];
    
    fprintf('..............................\n')
    fprintf('From the orbit below ISS to Debris %d \n', j+1);
    fprintf('From a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', a0/1e3, inc0*180/pi,...
        RAAN0*180/pi);
    fprintf('To a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
        Omega2*180/pi);
    fprintf('Delta v = %f km/s \n',  Tdv/1e3);
    fprintf('TOF = %f days \n', TOF/param.TU);
    
    %% Proximity operations at the target
    RAAN2 = Omega2 -param.k*(param.mu/af)^3.5*cos(incf)*param.RDV1;
    Total_TOF = Total_TOF +param.RDV1/param.TU;
   
    fprintf('..............................\n')
    fprintf(' Proximity operations at the Debris %d \n', j+1)
    fprintf('From a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
        Omega2*180/pi);
    fprintf('To a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
        RAAN2*180/pi);
    fprintf('TOF = %f days \n', param.RDV1/param.TU);


    %% Leg 3
    param.legno =  param.legno + 1; 
    
    param.m0 = mass_after +  m_debris(j+1);
    % from debris orbit 2
    a0 = debris(j+1,1); inc0 = debris(j+1,2); RAAN0 = RAAN2;
    
    % to deorbit orbit
    af = param.target_altitude + param.Re; incf = inc0;
    
    [TOF,Tdv,lastOmega,semiMajor, inclination, RAAN,tSeg,~,~,~,mass] =  ...
    kluver(a0,  af, inc0,incf,RAAN0,param);
    
    Total_TOF = Total_TOF +TOF/param.TU;
    Total_dV = Total_dV + Tdv(end)/1e3;
    
    fprintf('..............................\n')
    fprintf('From debris %d to orbit below ISS \n', j+1);
    fprintf('From a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', a0/1e3, inc0*180/pi,...
        RAAN0*180/pi);
    fprintf('To a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
        lastOmega*180/pi);
    fprintf('Delta v = %f km/s \n', Tdv(end)/1e3);
    fprintf('TOF = %f days \n', TOF/param.TU);
    
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
    yyaxis right;
    p = plot(tSeg/param.TU,Tdv);
    plot_latex(p, 'time(days)', '$\Delta v$ (m/s)','', '' ,{});
    
    

    mass_after = param.m0/exp(Tdv(end)/param.Isp/param.g0);
    mf= [mf, param.m0  - mass_after];
    
end


% Proximity operations to handover the debris to the sheperd.
Total_TOF = Total_TOF +param.RDV2/param.TU;
Omega_0_new = lastOmega -param.k*(param.mu/af)^3.5*cos(incf)*param.RDV2;

fprintf('..............................\n')
fprintf('Proximity operations to handover the debris %d to the sheperd. \n', i)
fprintf('From a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
    lastOmega*180/pi);
fprintf('To a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
    Omega_0_new*180/pi);
fprintf('TOF = %f days \n', param.RDV2/param.TU);


fprintf('-------------------------------------- \n')
fprintf('Total Delta v = %f m/s \n',Total_dV);
fprintf('Total TOF = %f days \n', Total_TOF);
fprintf('Total fuel consumed = %f kg \n', sum(mf));
