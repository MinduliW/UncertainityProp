function [rr, vv, E, mass] = getPosition(time, fname, mu)

load('DebrisData.mat', 'DebrisData');

index = find({DebrisData(:).name} == string(fname));
 
mass = DebrisData(index).mass;


Norad1 = '33500';
Norad2 = '39766';
Norad3 = '33492';

SatID1 = num2str(-100000-str2double(Norad1));
SatID2 = num2str(-100000-str2double(Norad2));
SatID3 = num2str(-100000-str2double(Norad3));

if strcmp(fname,'H2AF15')
    SatID = SatID1;
elseif strcmp(fname,'ALOS2')
    SatID = SatID2;
elseif strcmp(fname,'GOSAT')
    SatID = SatID3;
end


newstr = strcat('JD', {' '}, num2str(time));

datestring = cspice_str2et(newstr);

cspice_et2utc(datestring,'C', 6);

State1 = cspice_spkezr(SatID, datestring, 'J2000', 'none', '399');

rr = State1(1:3)*1e3; 
vv = State1(4:6)*1e3; 

E = CoordConv.vec2orbElem(rr,vv,mu);

% [a,ecc,incl,RAAN,argp,nu] = ijk2keplerian(rr, vv);
% 
% E= [a, ecc, incl*pi/180, RAAN*pi/180, argp*pi/180, nu*pi/180];

end

% 
% function [rr, vv, E,epoch, mass] = getPosition(time, fname, mu, J2, req)
% 
% 
% % E includes true anomaly, not mean anomaly
% load('/Users/minduli/Astroscale_ADR/Main/common/DebrisData.mat', 'DebrisData');
% 
% % find the index of the corresponding debris
% index = find({DebrisData(:).name} == string(fname));
% 
% ephemData = DebrisData(index).Kepl;
% 
% epoch  = DebrisData(index).epoch;
% mass = DebrisData(index).mass;
% 
% E(1:5) = ephemData(1:5);
% 
% M0 = ephemData(6);
% 
% % compute theta
% n = sqrt(mu/E(1)^3);
% 
% M = wrapTo2Pi(n*(time -epoch)*86400+M0);
% 
% % newton to compute E
% E0=M;
% 
% for k=1:1:10
%     EF=E0-(E0-E(2)*sin(E0)-M)/(1-E(2)*cos(E0));
%     E0=EF;
% end
% 
% % finally computing theta
% theta=2*atan(sqrt((1+E(2))/(1-E(2)))*tan(E0/2.0));
% E(6)=theta;
% 
% p = E(1)*(1-E(2)^2);
% E(4) = wrapTo2Pi(E(4) -3/2*J2*(req/p)^2*n*cos(E(3))*((time -epoch)*86400));
% E(5) = wrapTo2Pi(E(5) +3/4*J2*(req/p)^2*n*(5*(cos(E(3)))^2-1)*((time -epoch)*86400));
% 
% % computing cartesian coordinates
% 
% r  = E(1)*(1-E(2)^2)/(1+E(2)*cos(E(6)));
% rx = r*cos(E(6));
% ry = r*sin(E(6));
% 
% rrloc = [rx; ry; 0];
% 
% R(1,1) = cos(E(4))*cos(E(5))-sin(E(4))*sin(E(5))*cos(E(3));
% R(1,2) = -cos(E(4))*sin(E(5))-sin(E(4))*cos(E(5))*cos(E(3));
% R(1,3) = sin(E(4))*sin(E(3));
% R(2,1) = sin(E(4))*cos(E(5))+cos(E(4))*sin(E(5))*cos(E(3));
% R(2,2) = -sin(E(4))*sin(E(5))+cos(E(4))*cos(E(5))*cos(E(3));
% R(2,3) = -cos(E(4))*sin(E(3));
% R(3,1) = sin(E(5))*sin(E(3));
% R(3,2) = cos(E(5))*sin(E(3));
% R(3,3) = cos(E(3));
% 
% rr = R*rrloc;
% 
% thetap = sqrt(E(1)*(1-E(2)^2)*mu)/r^2;
% rp = E(1)*E(2)*(1-E(2)^2)*thetap*sin(E(6))/(1+E(2)*cos(E(6)))^2;
% 
% xp = rp*cos(E(6))-r*thetap*sin(E(6));
% yp = rp*sin(E(6))+r*thetap*cos(E(6));
% 
% vvloc = [xp;yp;0];
% vv = R*vvloc;
% rr = rr;
% E(1) = E(1);