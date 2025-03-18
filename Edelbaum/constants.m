%% constants 
param.TU = 86400;  % s
param.mu = 3.986005e14;%m^2/s^3
param.Re = 6378137; % m
param.J2 = 1.08262668e-3;
param.N = 25;
param.k = 3*param.J2*param.Re^2/(2*param.mu^3); %
param.AU =  149597870691; %m
param.Rs = 695800e3; %m
param.g0 = 9.80665; %m/s


% % fsolve options
% options = optimoptions('fsolve','Display','iter');
% options.Algorithm = 'levenberg-marquardt';
