% just calculating the nonlinear indicies. 
clear,clc;

ppLEO;
periodInit = 2*pi*sqrt(a0^3/param.mu);
param.tf = 15*periodInit;

[rr, vv] = CoordConv.po2pv([a0,  0.0064, deg2rad(99.2248),  deg2rad(5), ...
    deg2rad(5), deg2rad(5)],param.mu);

% variance 
r_var = 1e3/param.LU ;
v_var = 1/param.LU*param.TU; 

CartX = [rr;vv];
sigmaCart = [r_var, r_var, r_var, v_var, v_var, v_var];

mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/CovarianceNLI.cpp -ldace 
params = [param.mu ,param.J2, param.Re, param.tf];

[nlIndicies]= CovarianceNLI(params, CartX , sigmaCart');
