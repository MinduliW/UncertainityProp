function [Taeio,etaavg] = thrustGuidance(vv, a,e,i,Om,aop, nu, param)


E = 2*atan(tan(nu/2*sqrt(1-e)/sqrt(1+e)));
    
arglat = nu + aop;

%% Eccentricity
alphae = atan2(sin(nu),(cos(nu)+cos(E)));
etae = 1/2*(1+2*e*cos(nu) + cos(nu)^2)/(1+e*cos(nu));
betae = 0;


Te = (etae>param.etae)*sign(0-e)...
    *[cos(betae)*sin(alphae);cos(betae)*cos(alphae);sin(betae)];

if (isnan(norm(Te)))
    Te = [0;0;0];
    
end


%% Semi major axis 
alphaa = atan(e*sin(nu)/(1+e*cos(nu)));
betaa = 0;
etaa = norm(vv)*sqrt(a/param.mu*(1-e)/(1+e));

Ta = (etaa>param.etaa)*sign(param.smaexp-a)*[cos(betaa)*sin(alphaa);cos(betaa)*cos(alphaa);sin(betaa)];

if isnan(norm(Ta))
    Ta = [0;0;0];
end


%% Inclination

alphai = 0;
betai = sign(cos(arglat))*pi/2;
etai = abs(cos(arglat))/(1+e*cos(nu))*(sqrt(1-e^2*sin(aop)^2)-e*abs(cos(aop)));
Ti = (etai>param.etai)*sign(param.incexp-i)*[cos(betai)*sin(alphai);cos(betai)*cos(alphai);sin(betai)];

if (isnan(norm(Ti)))
    Ti = [0;0;0];
end

%% RAAN
alphaOm = 0; 
betaOm = sign(sin(arglat))*pi/2;

dRAAN = acos(cos(Om - param.RAANexp));

% if distance < 1e-3
%     distance = 0;
% end
sgndistance = -sin(Om - param.RAANexp)/(1 - cos(Om - param.RAANexp)^2)^(1/2);

etaRAAN = abs(sin(arglat))/(1+e*cos(nu))*(sqrt(1-e^2*cos(aop)^2)-e*abs(sin(aop)));

TOm = (etaRAAN>param.etaRAAN)*...
    sign(sgndistance)*[cos(betaOm)*sin(alphaOm);cos(betaOm)*cos(alphaOm);sin(betaOm)];


if (isnan(norm(TOm)))
    TOm = [0;0;0];
end

% 

%% AOP


% dAOP = acos(cos(aop - param.AOPexp));
% alphaAOP = atan((1+e*cos(nu))/(2+ e*cos(nu))*cot(nu));
% denom = sin(alphaAOP-nu)*(1+ e*cos(nu))-cos(alphaAOP)*sin(nu);
% betaAOP = atan(e*cot(i)*sin(arglat)/denom);
% etaAOP = abs(cos(arglat))/(1+e*cos(nu))*(sqrt(1-e^2*sin(aop)^2)-e*abs(cos(aop)));
% TAOP = (etaAOP>param.etaAOP)*sign(param.AOPexp-i)*...
%     [cos(betaAOP)*sin(alphaAOP);cos(betaAOP)*cos(alphaAOP);sin(betaAOP)];
% 
% if (isnan(norm(TAOP)))
%     TAOP = [0;0;0];
% end
% 



%%


Taeio = abs(param.smaexp-a)/100*Ta ...
    + abs(e)/1e-1*Te ...
    + abs(param.incexp-i)/(0.0005*pi/180)*Ti+ abs(dRAAN)/(5*pi/180)*TOm;




% if distance*180/pi > 1
%    distance*180/pi
%    param.RAANexp*180/pi
%    Om*180/pi
%     
% disp('--------------')
% 
% end

etaavg = (etae+etaa+etai+etaRAAN)/4;



if norm(Taeio)<1e-6  %|| etaavg<0.99       %param.dutyRatio
    Taeio = [0;0;0];
else
    Taeio = Taeio/norm(Taeio);
end

end