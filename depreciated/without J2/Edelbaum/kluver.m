function [TOF,deltaV, lastOmega, semiMajor, inclination, RAAN,tSeg,decayrate,beta,fs,mass] = ...
    kluver(a0, af, inc0,incf,Raan0, param)
% FUNCTION NAME:kluver(V0, Vf, inc0,incf,Raan0,param)
%
% DESCRIPTION:
%   Computes the TOF, total deltav and the destination RAAN according to
%   Craig A. Kluver's Edelbaum based method. Takes eclipses into account.
%
% INPUT: V0 (m/s), Vf (m/s),  inc0(rad), incf(rad), Raan0(rad)
%
% OUTPUT: TOF (s) , Total deltav (m/s), lastOmega (rad)
%
% Additional: mu (m^3/s^2), J2, Isp (s), g0 (m/s)

[semiMajor,inclination, ~, tSeg, deltaV_Edel,~,~,beta,f] =...
    Kechichan_Algorithm(a0, af, inc0,incf,param);

deltaV = 0.9766*deltaV_Edel;

%dVold = deltaV;
% Calculate velocities 
V0 = sqrt(param.mu / a0); %m/s
Vf = sqrt(param.mu/ af);%m/s


RAAN = zeros(size(tSeg));

if length(Raan0) >1
    Raan0 = Raan0(1);
end

RAAN(1) = Raan0;
decayrate(1) = 0;

fs  =[f]; 

% Go through each time segment
for i = 1:length(tSeg)-1
   
    if param.eclipses == true
        weight = shadowArc(semiMajor(i), inclination(i), RAAN(i), tSeg(i),param);
        weight = max([param.dutyRatio,weight]);
    else
        weight = param.dutyRatio;
    end
    
    
    % calculate the spacecraft mass vector.
    mk = param.m0/exp(deltaV(i)/param.Isp/param.g0);
    
    %mk = param.m0 - (weight*param.T)/param.Isp/param.g0*(tSeg(i))

    mass(i) = mk;
    % calculate average thrust.
    fk = param.Tmax/mk;
    
    
    fs = [fs, fk];
    if (weight*fk) == 0
        tSeg(i + 1)  =  tSeg(i);
    else
        
        tSeg(i + 1)  =  tSeg(i) + (deltaV(i+1) - deltaV(i))/(weight*fk);
    end
    
  
    if param.drag == true
        % Drag acceleration at tk
        f_d_k = drag_acceleration(semiMajor(i), param.t0+tSeg(i)/param.TU, param);
        % sma
        a_d_k  = param.mu./(V0^2 + f_d_k^2*tSeg(i).^2 -2*V0*f_d_k*tSeg(i)); %m
        
        % drag acceleratoin at t_k+1
        f_d_kplus1 =drag_acceleration(semiMajor(i+1), param.t0 + tSeg(i+1)/param.TU, param);
        a_d_kplus1  = param.mu./(V0^2 + f_d_kplus1^2*tSeg(i+1).^2 -2*V0*f_d_kplus1*tSeg(i+1)); %m
        
        % Calculate the percentage change of a. 
        a_drag_pchange = abs(a_d_k - a_d_kplus1)/a_d_kplus1;
        
        decayrate(i+1) = 0;
       if a_drag_pchange > param.a_pchangeLimit && abs(a_d_k - a_d_kplus1)<semiMajor(i)
            
            semiMajorold =  semiMajor(i+1);
            semiMajor(i+1) = semiMajor(i+1)-abs(a_d_k - a_d_kplus1);
            decayrate(i+1) = (semiMajorold - semiMajor(i+1))/(tSeg(i+1) - tSeg(i));
            
            % Redoing Edelbaum from i+1 to end. 
            Deltai = abs(incf - inclination(i+1));
            Vstart = sqrt(param.mu / semiMajor(i+1));
            DV_tospend= 0.9766*sqrt(Vstart^2 + Vf^2 - 2*Vstart*Vf*...
                cos(pi/2*Deltai));
            TOFtospend = abs(DV_tospend)/fk; %s
            future_time = linspace(0,TOFtospend, param.N -i); %s
            deltaV_newarray = future_time*fk;
            tSeg(i+1:end)=  tSeg(i+1) + future_time;
            deltaV(i+1:end) = deltaV(i+1) + deltaV_newarray;
      
      end
    end
        
%     avgsma = (semiMajor(i)+ semiMajor(i+1))/2; 
%     avginc =  (inclination(i)+ inclination(i+1))/2; 
%     
    n = sqrt(param.mu./semiMajor(i).^3);  %1/s
    
    RAANdot = -3/2*param.J2*n.*(param.Re./semiMajor(i)).^2.*cos(inclination(i));
    
    RAAN(i+1) =  RAAN(i) + RAANdot.*(tSeg(i+1) - tSeg(i));
    
    
    if ~isreal(RAAN(i+1))
        RAANdot
        
    end
    
    
    
    
end


mk = param.m0/exp(deltaV(end)/param.Isp/param.g0);
      
mass(i+1) = mk;


TOF = tSeg(end);
Tdv = deltaV(end);
lastOmega = RAAN(end);
end
