function [dx,coemean,coe] = forwardpropNothrust(t,x, tvals,param)

% Get A, B, C ,D and dB/dx for x.
[A,B,~, ~, ~] = getAandBMatrices(x,param);


% get velocity vector
coe = CoordConv.mee2coe(x);

[rr, vv] = CoordConv.po2pv(coe, param.mu);


%% Guidance adjustments 

% convert to hill
hill = kep2hill(coe, param.mu);

% calculate mean hill
hillMean = osculating2meanHill(hill, param.mu, param.J2, param.Re);

coemean = hill2kep(hillMean, param.mu);

[~, vvmean] = CoordConv.po2pv(coemean, param.mu);

if ~isreal(coemean)
    disp('Error: mean orbital elements are not real');
    
end
% 
% if param.eccAdj == true
%     
%     % Calculate the expected sma, inc and RAAN at this stage.
%     param.smaexp = interp1(tvals, param.Edelsma,t);
%     param.incexp = interp1(tvals, param.Edelinc,t);
%     param.RAANexp = interp1(tvals, param.EdelRAAN,t);
%     
%     
%     Emean = 2*atan(tan(coemean(6)/2*sqrt(1-coemean(2))/sqrt(1+coemean(2))));
%     
%     unitGuidance = thrustGuidance(vvmean,coemean(1),coemean(2),coemean(3),...
%         coemean(4),coemean(6),Emean,coemean(5), param);
%    
%     unitEdel = [0; 0; 0];
%     
%     if norm(unitEdel + unitGuidance) ~= 0
%         % get the unit combined thrust
%         unitvector = (1-param.eps)*unitEdel + param.eps*unitGuidance;
%         unitTotalthrust = unitvector/norm(unitvector);
%         Thrust = f*eta*unitTotalthrust;
%         
%         
%     end
%     
% end



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
    drag_term = 0;
end

 
%% State derivates

dx(1:6,1) = A+ B*deltaJ2;
dx(7,1) = 0;

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
C3 = 1-h^2 -k^2;
deltaJ2(1) =  -3/2*muJ2rE2r4*(1-12*C1^2/s4);
deltaJ2(2) =  -12*muJ2rE2r4*C1*C2/s4;
deltaJ2(3) =  -6*muJ2rE2r4*(C3*C1/s4);


