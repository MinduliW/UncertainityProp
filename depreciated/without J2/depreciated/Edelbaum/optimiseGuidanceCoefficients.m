function [difference,x] = optimiseGuidanceCoefficients(x,x0, betaedel,f, t2, param)


param.eps = x(1); 
param.etae = x(2);
param.etaa = x(3);
param.etai = x(4);

options = odeset('RelTol',1e-5);

period = 2*pi*sqrt(x0(1)^3/param.mu);

[time,states] = ode45(@(t,x) forwardPropagator(t,x, betaedel,f, t2, param),...
    [0,t2(end)],[x0';param.m0],options);


smadiff = zeros(size(time));
incdiff = zeros(size(time));
Omegadiff = zeros(size(time));

eccdiff = zeros(size(time));
for i = 1:length(time)
    [~, ~, ~, ~, ~, coemean,~] =...
        forwardPropagator(time(i)/param.TU,states(i,1:7)',...
        betaedel,f , t2, param);
    
    % get the sma, inclination and raan from edelbaum
    exptsma = interp1(t2, param.Edelsma, time(i));
    exptinc = interp1(t2, param.Edelinc, time(i));
    exptRAAN = interp1(t2, param.EdelRAAN, time(i));
    
    smadiff(i)= abs(coemean(1)-exptsma); % /coords(end,1)*100;
    incdiff(i) = abs(coemean(3)-exptinc)*180/pi;
    Omegadiff(i) = abs(wrapToPi(mod(coemean(4)- exptRAAN, 2*pi)))*180/pi;
    eccdiff(i) = abs(coemean(2));
    
    
end



difference = sum(smadiff)/mean(exptsma)+100*sum(incdiff)+sum(Omegadiff)+ sum(eccdiff);

end
