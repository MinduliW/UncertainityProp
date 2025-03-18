function [c,ceq] = nonlinear_constraint_general(x, debris,m_debris, param)

totalSmaChange = 0;
totalIncChange = 0;
%param.t0 = x(end); 
datevec = '19-June-2024';
param.t0 = juliandate(datevec,'dd-mmm-yyyy'); %+ 1000;

% Leg 1 : From debris orbit 1 to 350 km orbit, No RAAN matching required. 
param.m0  = param.m_dry_servicer + param.m_propellant + m_debris(1);

% from debris orbit 
a0 = debris(1,1);inc0 = debris(1,2);RAAN0 = debris(1,3); 

% to deorbit orbit
af = param.target_altitude + param.Re; incf = inc0;

totalSmaChange = totalSmaChange + abs(af - a0);
totalIncChange = totalIncChange + abs(incf -inc0);


[TOF,Tdv, lastOmega,~, ~, ~,~]=kluver(a0,  af, inc0,incf,RAAN0,param);

Total_TOF = TOF/param.TU; Total_dV = Tdv/1e3; 

mass_after = param.m0/exp(Tdv/param.Isp/param.g0); 

wt = zeros(1,length(m_debris)-1);

k = 1; 
currstart = 1; 
for j = 1: length(m_debris)-1
    %% Proximity operations to handover the debris to the sheperd.
    Omega_0_new = lastOmega -param.k*(param.mu/af)^3.5*cos(incf)*param.RDV2;
    
    Total_TOF = Total_TOF + param.RDV2/param.TU;
   
    %% From deorbit orbit to debris orbit, RAAN matching required.
    
    param.m0 = mass_after -  m_debris(j);
    
    % from deorbit orbit
    a0 = param.target_altitude + param.Re; inc0 = incf; RAAN0 = Omega_0_new;
    
    % to debris orbit 2
    af =  debris(j+1,1);
    incf = debris(j+1,2);
    RAANf_t0 = debris(j+1,3)-param.k*(param.mu/af)^3.5*cos(incf)*Total_TOF*param.TU;
    if param.optimiseRAAN == true
      RAANf_t0 = (x(3*j))-param.k*(param.mu/af)^3.5*cos(incf)*Total_TOF*param.TU;
    end
    
    
        
    totalSmaChange = totalSmaChange + abs(af - a0);
    totalIncChange = totalIncChange + abs(incf -inc0);

    
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


    %% Proximity operations at the target
    RAAN2 = Omega2 -param.k*(param.mu/af)^3.5*cos(incf)*param.RDV1;
    Total_TOF = Total_TOF +param.RDV1/param.TU;
   
    % Leg 3 
    param.m0 = mass_after +  m_debris(j+1);
    % from debris orbit 2
    a0 = debris(j+1,1); inc0 = debris(j+1,2); RAAN0 = RAAN2;
    
    % to deorbit orbit
    af = param.target_altitude + param.Re; incf = inc0;
    
    totalSmaChange = totalSmaChange + abs(af - a0);
    totalIncChange = totalIncChange + abs(incf -inc0);

    
    [TOF,Tdv,lastOmega] =  kluver(a0,  af, inc0,incf,RAAN0,param);
    
    Total_TOF = Total_TOF +TOF/param.TU;
    Total_dV = Total_dV + Tdv/1e3; 
  
    mass_after = param.m0/exp(Tdv/param.Isp/param.g0);
    

end


% Proximity operations to handover the debris to the sheperd.
Total_TOF = Total_TOF +param.RDV2/param.TU;
    
  
c = [];
ceq = [Total_TOF - 1e3];
ceq =[];
