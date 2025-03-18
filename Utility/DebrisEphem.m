function [rr, vv, E] = DebrisEphem(epoch, num, ephemData, mu, J2, req)

num = num+1; % for consistency with ESA

if nargin<3
    ephemData = load('Debris.txt');
end

E(1:5) = ephemData(num, 3:7);
E(1) = E(1)/1000;

M0 = ephemData(num,8);
dateM = ephemData(num, 2);

% compute theta
n = sqrt(mu/E(1)^3);
M = wrapTo2Pi(n*(epoch-dateM)*86400+M0);

% newton to compute E
E0=M;

for k=1:1:10
    EF=E0-(E0-E(2)*sin(E0)-M)/(1-E(2)*cos(E0));
    E0=EF;
end

% finally computing theta
theta=2*atan(sqrt((1+E(2))/(1-E(2)))*tan(E0/2.0));
E(6)=theta;

%% CHANGE MADE HERE
E(6) = M0;
%%


p = E(1)*(1-E(2)^2);
E(4) = wrapTo2Pi(E(4) -3/2*J2*(req/p)^2*n*cos(E(3))*((epoch-dateM)*86400));
E(5) = wrapTo2Pi(E(5) +3/4*J2*(req/p)^2*n*(5*(cos(E(3)))^2-1)*((epoch-dateM)*86400));

% computing cartesian coordinates

r  = E(1)*(1-E(2)^2)/(1+E(2)*cos(E(6)));
rx = r*cos(E(6));
ry = r*sin(E(6));

rrloc = [rx; ry; 0];

R(1,1) = cos(E(4))*cos(E(5))-sin(E(4))*sin(E(5))*cos(E(3));
R(1,2) = -cos(E(4))*sin(E(5))-sin(E(4))*cos(E(5))*cos(E(3));
R(1,3) = sin(E(4))*sin(E(3));
R(2,1) = sin(E(4))*cos(E(5))+cos(E(4))*sin(E(5))*cos(E(3));
R(2,2) = -sin(E(4))*sin(E(5))+cos(E(4))*cos(E(5))*cos(E(3));
R(2,3) = -cos(E(4))*sin(E(3));
R(3,1) = sin(E(5))*sin(E(3));
R(3,2) = cos(E(5))*sin(E(3));
R(3,3) = cos(E(3));

rr = R*rrloc;

thetap = sqrt(E(1)*(1-E(2)^2)*mu)/r^2;
rp = E(1)*E(2)*(1-E(2)^2)*thetap*sin(E(6))/(1+E(2)*cos(E(6)))^2;

xp = rp*cos(E(6))-r*thetap*sin(E(6));
yp = rp*sin(E(6))+r*thetap*cos(E(6));

vvloc = [xp;yp;0];
vv = R*vvloc;
rr = rr;
E(1) = E(1);