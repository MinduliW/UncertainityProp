function writeparams(param,N)

% FUNCTION NAME:
%   writeparams
%
% DESCRIPTION:
%   writes the SCVX problem parameters to a text file to be retrieved by
%   CPP
%
% INPUT:
%   param - (struct []) Problem parameter
%   N - (double) number of nodes 


tvec = linspace(param.t0, param.tf,N+1); 
dt = tvec(2)-tvec(1);

fileID = fopen('params.dat','w');
fprintf(fileID,'%20.20f \n',param.LU);
fprintf(fileID,'%20.20f \n',param.TU);
fprintf(fileID,'%20.20f \n',param.tf);
fprintf(fileID,'%20.20f \n',param.Tmax);
fprintf(fileID,'%20.20f \n',N);
fprintf(fileID,'%20.20f \n',param.m0);
fprintf(fileID,'%20.20f \n',param.mu);
fprintf(fileID,'%20.20f \n',param.Re);
fprintf(fileID,'%20.20f \n',dt);
fprintf(fileID,'%20.20f \n',param.method);
fprintf(fileID,'%20.20f \n',param.Isp);
fprintf(fileID,'%20.20f \n',param.g0);
fprintf(fileID,'%20.20f \n',1.0);
fprintf(fileID,'%20.20f \n',param.Rs);
fprintf(fileID,'%20.20f \n',param.eclipse);
fprintf(fileID,'%20.20f \n',param.nu);
fprintf(fileID,'%20.20f \n',param.Area);
fprintf(fileID,'%20.20f \n',param.Cd);
fprintf(fileID,'%20.20f \n',param.MU);
fclose(fileID);
end
