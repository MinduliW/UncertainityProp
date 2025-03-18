function [dx,eta,Thrustedel, Thrust, beta,coemean,coe,etaavg] = forwardPropagator...
    (t,x, betas,fs , tvals, param)


% say collision at 1 day mark.



% Get A, B, C ,D and dB/dx for x.
[A,B,~, ~, ~] = getAandBMatrices(x,param);


coe = CoordConv.mee2coe(x);

[rr, vv] = CoordConv.po2pv(coe, param.mu);
   
f =  param.T/x(7); %  interp1(tvals, fs,t); %

% convert to hill
hill = kep2hill(coe, param.mu);

% calculate mean hill
hillMean = osculating2meanHill(hill, param.mu, param.J2, param.Re);

coemean = hill2kep(hillMean, param.mu);

[~, vvmean] = CoordConv.po2pv(coemean, param.mu);

% Eclipses
if param.eclipses == true
    eta = isEclipse(t, coe(1),coe(2), coe(3), coe(4), coe(5), coe(6), param);
else
    eta = 1;
end
    
ArgofLat = wrapTo2Pi(coe(5) + coe(6));
    

if param.dutyRatio ~= 1
if (ArgofLat>pi/2 && ArgofLat<=3*pi/2)
    eta = 0;
    
else
    eta = 1;
end
end

    
if param.Thrust == true
    %% Thrust
    % get velocity vector
    
    
    if  param.incincrease == true
        if (ArgofLat>pi/2 && ArgofLat<=3*pi/2)
            betas = -abs(betas);
            
        else
            betas = abs(betas);
        end
    else
        if (ArgofLat>=pi/2 && ArgofLat<3*pi/2)
            betas = abs(betas);
        else
            betas = -abs(betas);
        end
        
    end
    
    beta = interp1(tvals, betas,t);

    Thrust = eta*f*[0; cos(beta); sin(beta)];
    Thrustedel = Thrust;
else
    Thrustedel =[0;0;0];
    Thrust = [0;0;0];
    beta = 0;
end



%% Guidance adjustments 



if ~isreal(coemean)
    disp('Error: mean orbital elements are not real');
    
end

if param.guidance == true
    
    % Calculate the expected sma, inc and RAAN at this stage.
    param.smaexp = interp1(tvals, param.Edelsma,t);
    param.incexp = interp1(tvals, param.Edelinc,t);
    param.RAANexp = interp1(tvals, param.EdelRAAN,t);

    [unitGuidance,etaavg] = thrustGuidance(vvmean,coemean(1),coemean(2),coemean(3),...
        coemean(4),coemean(5), coemean(6), param);
    
    if param.Thrust == true
        unitEdel = [0; cos(beta); sin(beta)];
    else
        unitEdel = [0; 0; 0];
    end
    
    
    if norm(unitEdel + unitGuidance) ~= 0
        % get the unit combined thrust
        unitvector = (1-param.eps)*unitEdel + param.eps*unitGuidance;
        if norm(unitvector) == 0
            Thrust = [0;0;0];
        else
            unitTotalthrust = unitvector/norm(unitvector);
            Thrust = f*eta*unitTotalthrust;
            
        end
        
    end
    
    
%     if  etaavg <  param.etaavg%76 %556
%         Thrust = [0;0;0];
%     end
%     
    
end




%% J2
deltaJ2= J2acceleration(x,coe, param)';

%% Drag

if param.drag == true
    rho = param.dragfactor*Density_HP(coe(1) - param.Re);
    
    v = sqrt(param.mu/coe(1));
    a_d = 0.5*rho*param.Cd*param.Area*v^2/x(7);
    
    % Drag acceleration is in the opposite direction to the velocity.
    drag_accel_ECI = -a_d*vv/norm(vv);
    
    % obtain the ECI to RTN matrix.
    hhat = cross(rr, vv)/norm(cross(rr, vv));
    rhat = rr/norm(rr);
    
    % Now convert drag_accel_ECI to drag_accel_RTN
    drag_accel_RTN = [rhat'; cross(hhat, rhat)' ; hhat']*drag_accel_ECI;
    
    drag_term = B*drag_accel_RTN;
else
    drag_term = [0;0;0;0;0;0];
end





% tlimit = param.TU - 10*60;
% % switch off thrust for 30 mins. 
% if t > tlimit && t < tlimit+30*60
%     Thrust = [0;0;0];
% end

thrust_term = B*Thrust;

%% State derivates
%thrust_term

dx(1:6,1) = A +drag_term + thrust_term+ B*deltaJ2;

if param.Thrust == true && sum(Thrust)~= 0
    dx(7,1) = - param.T*eta/(param.Isp*param.g0);
else
    dx(7,1) = 0;
end

function deltaJ2 = J2acceleration(x,coe, param)
% p = x(1);
% f = x(2);g = x(3);
h = x(4);k = x(5);L = x(6);
J2 = 1.08262668e-3;
r =  coe(1)*(1-coe(2))/(1 + coe(2)*cos(coe(6)));

cosL = cos(L);sinL = sin(L);

s2 = 1+h*h+k*k;
r2 = r^2;
muJ2rE2r4 = param.mu*J2*param.Re^2/r2^2;
s4 = s2^2;
C1 = h*sinL - k*cosL;
C2 = h*cosL + k*sinL;
C3 = 1-h^2-k^2;
deltaJ2(1) =  -3/2*muJ2rE2r4*(1-12*C1^2/s4);
deltaJ2(2) =  -12*muJ2rE2r4*C1*C2/s4;
deltaJ2(3) =  -6*muJ2rE2r4*(C3*C1/s4);
