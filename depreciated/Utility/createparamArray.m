function paramsArray = createparamArray(param,N)

paramsArray = [param.Tmax, param.m0, param.mu, param.J2, param.Re, param.method,...
   param.LU, param.Area, param.Cd, param.MU, N, param.tf , param.nu, param.drag ,param.coordSys];


end
