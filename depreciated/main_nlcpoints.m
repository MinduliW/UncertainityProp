%% nonlinear index comparison
clear,clc, close all;

% Note: this is phasing 
% TO DO: solve a standard convex problem using this
addpath('Edelbaum');
addpath('Dynamics');
addpath('Utility');
addpath('SGP4routines_NAIF');

addmosek;
cspice_furnsh('SGP4routines_NAIF/kernelmac.txt')


problemParam; 
param.Topt = true;


[driftorbit, Total_dv1,TOF,Omegadiff ,leg1,leg2, waitTimevec, targetRAAN, wait,~]...
    =RAANmatchingmethod(v0,vf, inc0, incf, param);


%% Generate Geqoe and Eqoe point clouds. 
% generate a cloud of points in Geqoe.
x0Kep = [leg1.a(1), 0, leg1.inc(1), wrapTo2Pi(leg1.RAAN(1)), 0,...
    wrapTo2Pi(sqrt(param.mu/leg1.a(1)^3)*leg1.t(1))];

% convert starting position to osculating 
x0KepOsc = mean2oscKep(x0Kep, param.mu, param.J2, param.Re);

[rr, vv] = CoordConv.po2pv(x0KepOsc,param.mu);
x0cartOsc=  [rr;vv];

Geqoex0 = RV2GEq(x0cartOsc(1:3),x0cartOsc(4:6),param.mu, param.J2, param.Re);

% generate a cloud of points. 
variance = abs(Geqoex0)*0.2;

GeqoeX = [];
for i = 1:100
    GeqoeX(i,:) = [Geqoex0 + randn(1,6).*variance,0,0,0];

end

% convert each coordinate to Eq. 
EqoeX = [];
for i = 1:100
    [R,V] = GEq2RV(GeqoeX(i,1:6), 1e-13,param.mu, param.J2, param.Re); 
    EqoeX(i,:) = [RV2CEq(R,V,param.mu, param.J2, param.Re), 0,0,0];
    
end

%% Calculate the nonlinearity indicies.


mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/nlc.cpp -ldace -lgsl -lfmt -lcblas


% orbital period initial
periodInit = 2*pi*sqrt(leg1.a(1)^3/param.mu);
% get the time duration of the optimization 
dtOpt = 1000*periodInit;
param.t0 = 0; 
param.tf = dtOpt;
param.nu = 1;
param.method = 2; 

param.coordSys = 1; 
paramArray = createparamArray(param,100);
[~, ~, ~,nlIndxGeq] = nlc(paramArray, GeqoeX);

param.coordSys = 2; 
paramArray = createparamArray(param,100);
[~, ~, ~,nlIndxEq] = nlc(paramArray, EqoeX);


figure; hold on;
ylabels = {'\nu', 'p_1', 'p_2', 'L', 'q_1', 'q_2'}; 
for i = 1:6
subplot(3,2,i); 
plot(nlIndxGeq(:,i),'b');hold on;
p = plot(nlIndxEq(:,i));
plot_latex(p, 'Node', ylabels(i), '', '', {'GEqOE', 'EqOE'});
end

