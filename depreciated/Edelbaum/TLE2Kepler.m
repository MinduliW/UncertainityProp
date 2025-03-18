function [name,OE]= TLE2Kepler(file)

constants;
% Open the TLE file and read TLE elements
fid = fopen(file, 'rb');
L2c = fscanf(fid,'%71c%',1);
L3c = fscanf(fid,'%70c%',1);
name = fscanf(fid,'%75c%');
% fprintf(L2c);
% fprintf([L3c,'\n']);
fclose(fid);
% Open the TLE file and read TLE elements
fid = fopen(file, 'rb');
L2 = fscanf(fid,'%d%6d%*c%5d%*3c%*2f%f%f%5d%*c%*d%5d%*c%*d%d%5d',[1,9]);
L3 = fscanf(fid,'%d%6d%f%f%f%f%f%f%f',[1,8]);
fclose(fid);

epoch = L2(1,4)*24*3600;        % Epoch Date and Julian Date Fraction
%Db    = L2(1,5);                % Ballistic Coefficient
inc   = L3(1,3)*pi/180;                % Inclination [deg]
RAAN  = L3(1,4)*pi/180;                % Right Ascension of the Ascending Node [deg]
e     = L3(1,5)/1e7;            % Eccentricity
w     = L3(1,6)*pi/180;                % Argument of periapsis [deg]
M     = L3(1,7)*pi/180;                % Mean anomaly [deg]
n     = L3(1,8);                % Mean motion [Revs per day]

% Orbital elements

a = (mu/(n*2*pi/(24*3600))^2)^(1/3);     % Semi-major axis [km]


% Six orbital elements
OE = [a e inc RAAN w M epoch];
% fprintf('\n a [km]   e      inc [rad]  RAAN [rad]  w[rad]    MA [rad] \n ')
% fprintf('%4.2f  %4.4f   %4.4f       %4.4f     %4.4f    %4.4f', OE);
% 
% 

end
