function [SunVec_AU,SunVec_LU] = getSunPosVec(JD,param)
% FUNCTION NAME:getSunPosVec(JD) 
%
% DESCRIPTION:
%   Compute the Sun position vector (Earth to Sun Vector)
%
% INPUT: (double) Julian Date
%
% OUTPUT: ([1,3]) Earth to Sun vector in AU and km
%

% Julian centuries from epoch
T_UT1 = (JD - 2451545.0)/36525;

% Mean longitude of the sun 
lambdaMSun = 280.460 + 36000.771*T_UT1;

MSun = 357.5291092 + 35999.05034*T_UT1;

lambdaEcliptic = lambdaMSun + 1.914666471*sind(MSun) + 0.019994643*sind(2*MSun);

r_Sun = 1.000140612 - 0.016708617*cosd(MSun) - 0.000139589*cosd(2*MSun);

epsilon = 23.439291 - 0.0130042*T_UT1;

SunVec_AU = [r_Sun*cosd(lambdaEcliptic); 
    r_Sun*cosd(epsilon)*sind(lambdaEcliptic);
    r_Sun*sind(epsilon)*sind(lambdaEcliptic)];



SunVec_km = 1.495978707e8*SunVec_AU;

SunVec_LU = SunVec_km*1e3/param.LU;
end

