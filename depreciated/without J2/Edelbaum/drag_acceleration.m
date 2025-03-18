function a_d = drag_acceleration(semiMajor, timeinJD, param)

% Get rho from height.
%rho = density_calc(semiMajor - param.Re);

%rho = param.dragfactor *Density_HP(semiMajor - param.Re);

rho = param.dragfactor *Density_HP(semiMajor*param.LU - param.Re*param.LU)*param.LU^3/param.MU;


%rho = Density_nrlmsise00(semiMajor - param.Re, timeinJD, param);

v = sqrt(param.mu/semiMajor);
% calculate ad.
%a_d = rho/(0.1570/param.Re)*Bstar*v^2;
a_d = 0.5*rho*param.Cd*param.Area*v^2/param.m0;

