function [RAANS,result,debrisData, propmass]  = estimOptRAAN(debris,m_debris, param)

%% Leg 1 
propmass = 0;

datevec = '19-June-2024';
param.t0 = juliandate(datevec,'dd-mmm-yyyy'); %+ 1000;

% Leg 1 : From debris orbit 1 to 350 km orbit, No RAAN matching required. 
param.m0  = param.m_wet_servicer + m_debris(1);

% from debris orbit 
a0 = debris(1,1);inc0 = debris(1,2);RAAN0 = debris(1,3); 

% to deorbit orbit
af = param.target_altitude + param.Re; incf = inc0;

[TOF,Tdv, lastOmega,~, ~, ~,~]=kluver(a0,  af, inc0,incf,RAAN0,param);

Total_TOF = TOF/param.TU; Total_dV = Tdv(end)/1e3; 

mass_after = param.m0/exp(Tdv(end)/param.Isp/param.g0); 

propmass = propmass+(param.m0- mass_after);

debrisData = []; 
k = 1; 
RAANS =[];
for j = 1: length(m_debris)-1
    %% Proximity operations to handover the debris to the sheperd.
    Omega_0_new = lastOmega -param.k*(param.mu/af)^3.5*cos(incf)*param.RDV2;
    
    Total_TOF = Total_TOF + param.RDV2/param.TU;
   
    %% From deorbit orbit to debris orbit, RAAN matching required.
    
    param.m0 = mass_after -  m_debris(j);
    
    % from deorbit orbit
    a0 = param.target_altitude + param.Re; inc0 = incf; RAAN0 = Omega_0_new;
    debrisData = [debrisData , sqrt(param.mu/a0)];
    debrisData = [debrisData , inc0];
    
    % to debris orbit 2
    af =  debris(j+1,1);
    incf = debris(j+1,2);
   
    [TOF,Tdv, Omega2,~, ~, ~,~]=kluver(a0,  af, inc0,incf,RAAN0,param);

    k = k +1; 
    
    Total_TOF = Total_TOF +TOF/param.TU;
    Total_dV = Total_dV + Tdv(end)/1e3; 
    
    % We need omega at the start of the whole trajectory.
    Omstt = Omega2 +param.k*(param.mu/af)^3.5*cos(incf)*Total_TOF*param.TU;
    RAANS = [RAANS, (Omstt)];
    
    debrisData = [debrisData ,  wrapTo2Pi(Omstt)];
    
    mass_after = param.m0 /exp(Tdv(end)/param.Isp/param.g0);

    propmass = propmass+(param.m0- mass_after);

    %% Proximity operations at the target
    RAAN2 = Omega2 -param.k*(param.mu/af)^3.5*cos(incf)*param.RDV1;
    Total_TOF = Total_TOF +param.RDV1/param.TU;
   
    % Leg 3 
    param.m0 = mass_after +  m_debris(j+1);
    % from debris orbit 2
    a0 = debris(j+1,1); inc0 = debris(j+1,2); RAAN0 = RAAN2;
    
    % to deorbit orbit
    af = param.target_altitude + param.Re; incf = inc0;
    
    [TOF,Tdv,lastOmega] =  kluver(a0,  af, inc0,incf,RAAN0,param);
    
    Total_TOF = Total_TOF +TOF/param.TU;
    Total_dV = Total_dV + Tdv(end)/1e3; 
  
    mass_after = param.m0/exp(Tdv(end)/param.Isp/param.g0);

    propmass = propmass+(param.m0- mass_after);

end
% Proximity operations to handover the debris to the sheperd.
Total_TOF = Total_TOF +param.RDV2/param.TU;
    
  
if param.Topt == true

  result = Total_TOF; 
else
   result = Total_dV;
end


 
end
