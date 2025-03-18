function [a,i, tf, t, deltav,DeltaVtotcheck,V_timehistory,beta,f] = ...
    Kechichan_Algorithm(a0, af, inc0,incf,param)

% Calculate Delta i 
Deltai = abs(incf- inc0); %rad 

% Constant thrust acceleration 
f = param.Tmax/ param.m0; % m/s^2

% Calculate velocities 
V0 = sqrt(param.mu / a0); %m/s
Vf = sqrt(param.mu/ af);%m/s

% Initial yaw steering angle 
beta0 = atan2(sin(pi/2*Deltai),(V0/Vf - cos(pi/2*Deltai)));

% compute total delta v 
%DeltaVtot = V0*cos(beta0) - V0*sin(beta0)/(tan(pi/2*Deltai + beta0)); %m/s

% Check.
if Deltai <= 2
    DeltaVtotcheck = sqrt(V0^2 + Vf^2 - 2*V0*Vf*cos(pi/2*Deltai)); %m/s
else
    DeltaVtotcheck = V0 + Vf;
end

% Transfer time with continous thrust. 
tf = abs(DeltaVtotcheck)/f; %s

% Obtaining time histories. 
t = linspace(0,tf, param.N); %s

% Time history of delta v 
deltav = f*t; %m/s 

% beta = SOLVE FOR BETA

beta = atan2(V0*sin(beta0), V0*cos(beta0)-f*t); 

a = param.mu./(V0^2 + f^2*t.^2 -2*V0*f*t*cos(beta0)); %m

% inclination
inside_atan = (f*t - V0*cos(beta0))/(V0*sin(beta0));

i = inc0+ sign(incf -inc0)*2/pi*(atan(inside_atan) + pi/2 -beta0); %rad

V_timehistory = sqrt(V0^2 - 2*V0*f*t*cos(beta0) +f^2*t.^2); %m/s

function y = atan3(a, b)

% four quadrant inverse tangent

% input

%  a = sine of angle
%  b = cosine of angle

% output

%  y = angle (radians; 0 =< c <= 2 * pi)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon = 0.0000000001;

pidiv2 = 0.5 * pi;

if (abs(a) < epsilon)
   y = (1 - sign(b)) * pidiv2;
   return;
else
   c = (2 - sign(a)) * pidiv2;
end

if (abs(b) < epsilon)
   y = c;
   return;
else
   y = c + sign(a) * sign(b) * (abs(atan(a / b)) - pidiv2);
end
