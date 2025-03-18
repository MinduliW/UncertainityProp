function [param] =SCvxParam(x0cartMEAN,xfcartMEAN, tvec, param) % 



param.Tmax = param.T; %30/param.MU/(param.LU)*param.TU^2;
% param.Isp = 285/param.TU; 


param.t0 = tvec(1);
param.tf = tvec(end); 

param.x0 = CoordConv.kepler2MEOE([x0cartMEAN(1), x0cartMEAN(2), x0cartMEAN(5),...
x0cartMEAN(4),x0cartMEAN(3),x0cartMEAN(6)]);

param.xf =CoordConv.kepler2MEOE([xfcartMEAN(1), xfcartMEAN(2), xfcartMEAN(5),...
xfcartMEAN(4),xfcartMEAN(3),xfcartMEAN(6)]);


% this is mean coordinate.
period1 = 2*pi*sqrt(param.x0(1)^3/param.mu);
period2 = 2*pi*sqrt(param.xf(1)^3/param.mu);



period = 0.5*(period1 + period2);
revs = round(param.tf/period);
param.xf(6) = wrapTo2Pi(param.xf(6))+ (revs-1)*2*pi;


% calculate number of nodes
% if revs > 108
%     N = round((param.tf-param.t0)/period*30);
% else
%     N = 200;
% end

%N = length(tvec);


param.tvec = tvec; %linspace(param.t0, param.tf,N+1); 
param.dt = param.tvec(2)-param.tvec(1);

fprintf("sma change (km) = %f \n" , abs(xfcartMEAN(1)-x0cartMEAN(1))*param.LU/1e3 );
fprintf("e change = %f \n" , abs(xfcartMEAN(2)-x0cartMEAN(2)) );
fprintf("inc change (deg)= %f \n" , abs(xfcartMEAN(3)-x0cartMEAN(3))*180/pi );
fprintf("RAAN change (deg)= %f \n" , wrapTo2Pi(xfcartMEAN(4)-x0cartMEAN(4))*180/pi  );
fprintf("AOP change = %f \n" , abs(xfcartMEAN(5)-x0cartMEAN(5))*180/pi  );

end
