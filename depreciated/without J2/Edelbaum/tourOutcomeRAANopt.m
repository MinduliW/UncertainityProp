function tourOutcomeRAANopt(x, debrisID, param)


noOfdrifts = length(debrisID) -1; 

RAANcoordinateInx = noOfdrifts*2+1; 

m_debris = zeros(size(debrisID));
debris = zeros(length(debrisID),6);

for i = 1:length(debrisID)
    [~, ~, debris(i,:), ~, m_debris(i)] = getPosition(0, ...
        debrisID(i),param.mu, param.J2, param.Re);
    
end


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

param.m0  = param.m_dry_servicer + param.m_propellant + m_debris(1);

% from debris orbit 
a0 = debris(1,1);inc0 = debris(1,3);RAAN0 = x(RAANcoordinateInx); 

% to deorbit orbit
af = param.target_altitude + param.Re; incf = inc0;

[TOF,Tdv, lastOmega,semiMajor, inclination, RAAN,tSeg]=kluver(a0,  af, inc0,incf,RAAN0,param);
mass_after = param.m0/exp(Tdv/param.Isp/param.g0); 


dVs(1) = Tdv/1e3;
TOFs(1) = TOF/param.TU;
mf(1) = param.m0  - mass_after;

fprintf('Leg 1 :From Debris 1 to orbit below ISS \n')
fprintf('From a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', a0/1e3, inc0*180/pi,...
    RAAN0*180/pi);
fprintf('To a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
    lastOmega*180/pi);
fprintf('Delta v = %f km/s \n', dVs(1));
fprintf('TOF = %f days \n', TOFs(1));

figure; 
subplot(3,1,1);hold on;
plot(tSeg/param.TU, a0*ones(size(tSeg)), '--');
plot(tSeg/param.TU, af*ones(size(tSeg)),'--');
p = plot(tSeg/param.TU,semiMajor);
plot_latex(p, 'time(days)', 'a (m)','', 'From Debris 1 to orbit below ISS' ,{'$a_0$','$a_f$','$a(t)$'});

subplot(3,1,2); hold on;
plot(tSeg/param.TU, rad2deg(inc0*ones(size(tSeg))), '--');
plot(tSeg/param.TU, rad2deg(incf*ones(size(tSeg))),'--');
p = plot(tSeg/param.TU,inclination*180/pi);
plot_latex(p, 'time(days)', 'inclination (deg)','', '' ,{'$i_0$','$i_f$','$i(t)$'});

subplot(3,1,3); hold on;
Omega_x0 = RAAN0 -param.k*(param.mu/a0)^3.5*cos(inc0)*tSeg;
plot(tSeg/param.TU, rad2deg(Omega_x0), '--');
p = plot(tSeg/param.TU,RAAN*180/pi);
plot_latex(p, 'time(days)', 'Omega (deg)','', '' ,{'$\Omega_0$','$\Omega(t)$'});


k = 1;
wt = [];
for i = 2: length(debrisID)

% Proximity operations to handover the debris to the sheperd.
mass_before = mass_after - m_debris(i -1);
Omega_0_new = lastOmega -param.k*(param.mu/af)^3.5*cos(incf)*param.RDV2;
TOFs = [TOFs , param.RDV2/param.TU];

fprintf('..............................\n')
fprintf('Proximity operations to handover the debris %d to the sheperd. \n', i-1)
fprintf('From a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
    lastOmega*180/pi);
fprintf('To a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
    Omega_0_new*180/pi);
fprintf('TOF = %f days \n', param.RDV2/param.TU);


%% From deorbit orbit to debris orbit 2, RAAN matching required. 
param.m0 = mass_before;

% from deorbit orbit
a0 = param.target_altitude + param.Re; inc0 = incf; RAAN0 = Omega_0_new;

% to debris orbit 2
af =  debris(i,1);incf = debris(i,3);
RAANf_t0 = x(RAANcoordinateInx+i-1)-param.k*(param.mu/af)^3.5*cos(incf)*(sum(TOFs)*param.TU);

endentry = i + (k -1);
startentry = endentry-1;

[Tdv,TOF,Omega,waittime] = propOpt(x(startentry:endentry), a0,inc0,RAAN0 , af, incf,...
    RAANf_t0, param);

wt(k) = waittime;
k = k+ 1;

dVs = [dVs, Tdv/1e3];
TOFs = [TOFs , TOF/param.TU];

mass_after = param.m0/exp(Tdv/param.Isp/param.g0); 
mf = [mf, param.m0- mass_after];

fprintf('..............................\n')
fprintf('From the orbit below ISS to Debris %d \n', i);
fprintf('From a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', a0/1e3, inc0*180/pi,...
    RAAN0*180/pi);
fprintf('To a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
    Omega*180/pi);
fprintf('Delta v = %f km/s \n', dVs(end));
fprintf('TOF = %f days \n', TOFs(end));



%% Proximity operations at the target
mass_before = mass_after + m_debris(i);
RAAN = Omega -param.k*(param.mu/af)^3.5*cos(incf)*param.RDV1;
TOFs = [TOFs , param.RDV1/param.TU];

fprintf('..............................\n')
fprintf(' Proximity operations at the Debris %d \n', i)
fprintf('From a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
    Omega*180/pi);
fprintf('To a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
    RAAN*180/pi);
fprintf('TOF = %f days \n', param.RDV1/param.TU);


%% Leg 3 

param.m0 = mass_before;

% from debris orbit 2
a0 = debris(i,1); inc0 = debris(i,3); RAAN0 = RAAN;

% to deorbit orbit
af = param.target_altitude + param.Re; incf = inc0;

[TOF,Tdv,lastOmega,semiMajor, inclination, RAAN,tSeg] =  kluver(a0,  af, inc0,incf,RAAN0,param);

dVs = [dVs, Tdv/1e3];
TOFs = [TOFs , TOF/param.TU];

mass_after = mass_before/exp(Tdv/param.Isp/param.g0);
mf = [mf, param.m0- mass_after];

fprintf('..............................\n')
fprintf('From debris %d to orbit below ISS \n', i);
fprintf('From a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', a0/1e3, inc0*180/pi,...
    RAAN0*180/pi);
fprintf('To a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
    lastOmega*180/pi);
fprintf('Delta v = %f km/s \n', dVs(end));
fprintf('TOF = %f days \n', TOFs(end));


figure; 
subplot(3,1,1);hold on;
plot(tSeg/param.TU, a0*ones(size(tSeg)), '--');
plot(tSeg/param.TU, af*ones(size(tSeg)),'--');
p = plot(tSeg/param.TU,semiMajor);
plot_latex(p, 'time(days)', 'a (m)','', 'From Debris to orbit below ISS' ,{'$a_0$','$a_f$','$a(t)$'});

subplot(3,1,2); hold on;
plot(tSeg/param.TU, rad2deg(inc0*ones(size(tSeg))), '--');
plot(tSeg/param.TU, rad2deg(incf*ones(size(tSeg))),'--');
p = plot(tSeg/param.TU,inclination*180/pi);
plot_latex(p, 'time(days)', 'inclination (deg)','', '' ,{'$i_0$','$i_f$','$i(t)$'});

subplot(3,1,3); hold on;
Omega_x0 = RAAN0 -param.k*(param.mu/a0)^3.5*cos(inc0)*tSeg;
plot(tSeg/param.TU, rad2deg(Omega_x0), '--');
p = plot(tSeg/param.TU,RAAN*180/pi);
plot_latex(p, 'time(days)', 'Omega (deg)','', '' ,{'$\Omega_0$','$\Omega(t)$'});



end


% Proximity operations to handover the debris to the sheperd.
TOFs = [TOFs , param.RDV2/param.TU];
Omega_0_new = lastOmega -param.k*(param.mu/af)^3.5*cos(incf)*param.RDV2;


fprintf('..............................\n')
fprintf('Proximity operations to handover the debris %d to the sheperd. \n', i)
fprintf('From a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
    lastOmega*180/pi);
fprintf('To a = %.2f km, i = %.2f deg and RAAN = %.2f deg \n', af/1e3, incf*180/pi,...
    Omega_0_new*180/pi);
fprintf('TOF = %f days \n', param.RDV2/param.TU);


dVT =  sum(dVs);
TOFT =  sum(TOFs);


fprintf('-------------------------------------- \n')
fprintf('Total Delta v = %f km/s \n',dVT);
fprintf('Total TOF = %f days \n', TOFT);
fprintf('Total fuel consumed = %f kg \n', sum(mf));

end