%% nonlinear index comparison
clear,clc; close all;

% Note: this is phasing 
% TO DO: solve a standard convex problem using this
addpath('Edelbaum');
addpath('Dynamics');
addpath('Utility');
addpath('SGP4routines_NAIF');

addmosek;
% cspice_furnsh('SGP4routines_NAIF/kernelmac.txt')

ppLEO;
%ppGeo; 

[rr, vv] = CoordConv.po2pv([a0,  0.0064, deg2rad(0.01),  deg2rad(5), ...
    deg2rad(5), deg2rad(5)],param.mu);

% variance 
r_var = 1e3/param.LU ;
v_var = 1/param.LU*param.TU; 

CartX(1,:) = [rr;vv];

% orbital period initial
periodInit = 2*pi*sqrt(a0^3/param.mu);
param.t0 = 0; 


dxCartMax = 3*[r_var, r_var,r_var, v_var , v_var, v_var];

MEEX0 =  CoordConv.vec2mee(CartX(1:3),CartX(4:6),param.mu);
GEqoeX0 = RV2GEq(CartX(1:3),CartX(4:6),param.mu, param.J2,param.Re);
CEqoeX0 = RV2CEq(CartX(1:3),CartX(4:6),param.mu, param.J2,param.Re);
kepx0   = vec2orbElem(CartX(1:3),CartX(4:6),param.mu);

% propagate for 15 orbits and get the absol. pos. 
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
%[t,y] = ode45(@(t,y) propagateCart(t,y, param), [0, param.tf], CartX, options);

% CartXend = y(end,:);
% 
% MEEXf   = CoordConv.vec2mee(CartXend(1:3),CartXend(4:6),param.mu);
% GEqoeXf = RV2GEq(CartXend(1:3),CartXend(4:6),param.mu, param.J2,param.Re);
% CEqoeXf = RV2CEq(CartXend(1:3),CartXend(4:6),param.mu, param.J2,param.Re);
% kepxf   = vec2orbElem(CartXend(1:3),CartXend(4:6),param.mu);
% 

mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/CovarianceProp.cpp -ldace

norbits = 2:10:100;
for i = 1: length(norbits)


    param.tf = norbits(i)*periodInit;
    paramArray = createparamArray(param);

    %valkep(i,:) = kepX - kepx0;
    val(i,:)  =  1/3*randn(1,6);

    [fxe, fx0 , nlc(i),fxeMEE, fx0MEEc, nlcMEE(i), fxeGeq, fx0Geqc, nlcGeq(i),...
       fxeCeq, fx0Ceqc, nlcCeq(i), fxeKep(i,:), fx0kepc, nlcKep(i)]= ...
        CovarianceProp(paramArray, CartX , val(i,:));

%     dCartXX(i,:) = fxe- fx0;
%     dMEEXX(i,:) = fxeMEE -fx0MEEc;
%     dGEqoeXX(i,:) = fxeGeq - fx0Geqc;
%     dCEqoeXX(i,:) = fxeCeq - fx0Ceqc;
%     dKepXX(i,:) = fxeKep(i,:)' - fx0kepc;
end


figure;
semilogy(norbits, nlc, '-o'); hold on; 
semilogy(norbits, nlcMEE, '-o');
semilogy(norbits, nlcGeq, '-o');
semilogy(norbits, nlcCeq, '-o');
p = semilogy(norbits, nlcKep, '-o');
ylim([1e-7,1e2])
plot_latex(p, 'No. of orbits', '\nu','', '' , ...
    {'Cartesian', 'MEE', 'GeqOE', 'CeqOE', 'COE'});

