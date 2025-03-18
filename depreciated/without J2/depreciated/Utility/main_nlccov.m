%% nonlinear index comparison
clear,clc, close all;

addpath('Utility');
addpath('SGP4routines_NAIF');

addmosek;
cspice_furnsh('SGP4routines_NAIF/kernelwin.txt')


problemParam; 
param.Topt = true;

nsamples = 100;

x0cartOsc= [ -0.547674680736674;0.912767983841398;0;0.132602340906191; 0.079563422475774; 0.959260135812878];

% variance 
r_var = 10e3/param.LU ;
v_var = 5/param.LU*param.TU; 

CartX = [];


CartX(1,:) =  x0cartOsc';
for i = 2:nsamples+1
    val = x0cartOsc + [randn(1,3),0,0,0]'.*[r_var*ones(3,1); zeros(3,1)] + [0,0,0,randn(1,3)]'.*[zeros(3,1); v_var*ones(3,1)];
    CartX(i,:) =  [val']; 
end


% convert each coordinate to Eq. 
GEqoeX = [];
for i = 1:nsamples+1
    GEqoeX(i,:) = [RV2GEq(CartX(i,1:3),CartX(i,4:6),param.mu, param.J2, param.Re)];
    
end

% convert each coordinate to Eq. 
CEqoeX = [];
for i = 1:nsamples+1
    CEqoeX(i,:) = [RV2CEq(CartX(i,1:3),CartX(i,4:6),param.mu, param.J2, param.Re)];
    
end



%% transfrom to Geqoe. 
mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/GequProp.cpp -ldace -lgsl -lfmt -lcblas
mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/CequProp.cpp -ldace -lgsl -lfmt -lcblas
mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/MEEProp.cpp -ldace -lgsl -lfmt -lcblas


% orbital period initial
periodInit = 2*pi*sqrt(a0^3/param.mu);
param.t0 = 0; 
param.tf = 10*periodInit;
param.nu = 1;
param.method = 2; 


param.coordSys = 1; 
paramArray = createparamArray(param,nsamples+1);

%% 

[~, ~,cartxp, Geqoex ] = GequProp(paramArray, CartX);
[~, ~,~, Ceqoex ] = CequProp(paramArray, CartX);

[~, ~,~, MEEx ] = MEEProp(paramArray, CartX);


figure; hold on; 

for i = 1:nsamples+1
    p = plot(MEEx(:,1)-MEEx(1,1), MEEx(:,6)- MEEx(1,6), 'r.');
end


GeqoexRV =[];
CeqoexRV = [];
% lets convert Geqoex and Ceqoex to cartesian
for i = 1: nsamples+1
    [R,V]= GEq2RV(Geqoex(i,:) ,1e-13,param.mu, param.J2, param.Re);
    GeqoexRV(i,:) = [R,V];

    [R,V]= CEq2RV(Ceqoex(i,:) ,1e-13,param.mu, param.J2, param.Re);
    CeqoexRV(i,:) = [R,V];
end



figure;  hold on; 
plot(CartX(:,1)-CartX(1,1), CartX(:,2)-CartX(1,2), 'k.');
for i = 1:nsamples+1
    p = plot(cartxp(:,1)-cartxp(1,1), cartxp(:,2)-cartxp(1,2), 'r.');
end

plot_latex(p, 'X[Re]', 'Y[re]','', '100 LEO orbits (Cartesian)' ,{'Initial uncertainity', 'Final uncertainity'});


figure; hold on;

plot(CartX(:,1)-CartX(1,1), CartX(:,2)-CartX(1,2), 'k.');
for i = 1:nsamples+1
    p = plot((CeqoexRV(:,1)-CeqoexRV(1,1)), (CeqoexRV(:,2)-CeqoexRV(1,2)), 'g.');
end
for i = 1:nsamples+1
    p = plot((GeqoexRV(:,1)-GeqoexRV(1,1)),(GeqoexRV(:,2)-GeqoexRV(1,2)), 'r.');
end

plot_latex(p, 'X[Re]', 'Y[re]','', '100 LEO orbits (Final Geqoe-red, Final Ceqoe-green, Initial - black)' ,{});

