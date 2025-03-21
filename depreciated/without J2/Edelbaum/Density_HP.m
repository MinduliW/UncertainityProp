%---------------------------------------------------------------------------
%
% Density_HP: Computes the atmospheric density for the modified
%             Harris-Priester model.
%
% Inputs:
%   Mjd_TT      Terrestrial Time (Modified Julian Date)
%               mean equator and equinox of J2000 (EME2000, ICRF)
%   r_tod       Satellite position vector in the inertial system [m]
% Output:
%   density     Density [kg/m^3]
%
%---------------------------------------------------------------------------
function density = Density_HP(height)



% Harris-Priester atmospheric density model parameters 
% Height [km], minimum density, maximum density [gm/km^3]

h = [...
   100.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0,...
   210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0,...
   320.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0,...
   520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0, 680.0, 700.0,...
   720.0, 740.0, 760.0, 780.0, 800.0, 840.0, 880.0, 920.0, 960.0, 1000.0];

c_min = [...
   4.974e+05, 2.490e+04, 8.377e+03, 3.899e+03, 2.122e+03, 1.263e+03,...
   8.008e+02, 5.283e+02, 3.617e+02, 2.557e+02, 1.839e+02, 1.341e+02,...
   9.949e+01, 7.488e+01, 5.709e+01, 4.403e+01, 3.430e+01, 2.697e+01,...
   2.139e+01, 1.708e+01, 1.099e+01, 7.214e+00, 4.824e+00, 3.274e+00,...
   2.249e+00, 1.558e+00, 1.091e+00, 7.701e-01, 5.474e-01, 3.916e-01,...
   2.819e-01, 2.042e-01, 1.488e-01, 1.092e-01, 8.070e-02, 6.012e-02,...
   4.519e-02, 3.430e-02, 2.632e-02, 2.043e-02, 1.607e-02, 1.281e-02,...
   1.036e-02, 8.496e-03, 7.069e-03, 4.680e-03, 3.200e-03, 2.210e-03,...
   1.560e-03, 1.150e-03                                            ];        

% Satellite height
height = height/1000;

density = interp1(h,c_min, height);

if height> h(end)
    density = c_min(end);
end


if height< h(1)
    density = c_min(1);
end

density = density * 1.0e-12;      % [kg/m^3]