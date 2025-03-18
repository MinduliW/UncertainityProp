%% nonlinear index comparison
clear,clc, close all;

% Note: this is phasing 
% TO DO: solve a standard convex problem using this
addpath('Edelbaum');
addpath('Dynamics');
addpath('Utility');
addpath('SGP4routines_NAIF');

addmosek;
cspice_furnsh('SGP4routines_NAIF/kernelwin.txt')


ppGeo; 
param.Topt = true;



x0cartOsc= [ -0.547674680736674;0.912767983841398;0;0.132602340906191; 0.079563422475774; 0.959260135812878];

% variance 
r_var = 1e3/param.LU ;
v_var = 0.1/param.LU*param.TU; 

CartCov = [r_var^2*eye(3), zeros(3,3); zeros(3,3), v_var^2*eye(3)];
CartX = [];

   CartX(1,:) =  [x0cartOsc']; 

for i = 2:101
    val = x0cartOsc + [randn(1,3),0,0,0]'.*[r_var*ones(3,1); zeros(3,1)] + [0,0,0,randn(1,3)]'.*[zeros(3,1); v_var*ones(3,1)];
    CartX(i,:) =  [val']; 
end


%% transfrom to Geqoe. 
mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/RV2Geq.cpp -ldace -lgsl -lfmt -lcblas


% orbital period initial
periodInit = 2*pi*sqrt(param.x0(1)^3/param.mu);
% get the time duration of the optimization 
dtOpt = 100*periodInit;
param.t0 = 0; 
param.tf = dtOpt;
param.nu = 1;
param.method = 2; 


param.coordSys = 1; 
paramArray = createparamArray(param,1);

[JRV2Geq] = RV2Geq(paramArray, CartX(1,:));

Covariance = JRV2Geq(1:6,1:6)*JRV2Geq(1:6,1:6)';











%convert to Geqoe. 
GeqoeX = [];
GeqoeX(1,:) = [RV2GEq(rr, vv,param.mu, param.J2, param.Re),0,0,0];

for i = 2:101
    GeqoeX(i,:) = [RV2GEq(CartX(i-1,1:3), CartX(i-1,4:6),param.mu, param.J2, param.Re),0,0,0];

end





% convert to eqoe
EqoeX = [];
EqoeX(1,:) = [RV2CEq(rr, vv,param.mu, param.J2, param.Re),0,0,0];

for i = 2:101
    EqoeX(i,:) = [RV2CEq(CartX(i-1,1:3), CartX(i-1,4:6),param.mu, param.J2, param.Re),0,0,0];

end


mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/RV2Geq.cpp -ldace -lgsl -lfmt -lcblas
mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/RV2Ceq.cpp -ldace -lgsl -lfmt -lcblas





% orbital period initial
periodInit = 2*pi*sqrt(leg1.a(1)^3/param.mu);
% get the time duration of the optimization 
dtOpt = 15*periodInit;
param.t0 = 0; 
param.tf = dtOpt;
param.nu = 1;
param.method = 2; 


param.coordSys = 1; 
paramArray = createparamArray(param,101);

[J_Geq2RV] = Geq2RV(paramArray, GeqoeX);
[J_Eq2RV] = Ceq2RV(paramArray, EqoeX);


%CartCov = [r_var^2*eye(3), zeros(3,3); zeros(3,3), v_var^2*eye(3)];

GeqCov = J_Geq2RV(1:6,1:6)\CartCov*J_Geq2RV(1:6,1:6);
CeqCov = J_Eq2RV(1:6,1:6)\CartCov*J_Eq2RV(1:6,1:6);

[As, ~, ~,nlIndxGeq] = nlc(paramArray, GeqoeX);

% propagate the covariances. 
for i = 1:101
    ubd = 6*i;
    lbd = ubd - 5;

    A = As(lbd:ubd,1:6);

    %GeqCovProp(:,:.i) = 

end

GeqCovProp = A*GeqCov/A; 


param.coordSys = 2; 
paramArray = createparamArray(param,101);
[As, ~, ~,nlIndxEq] = nlc(paramArray, GeqoeX);

% propagate the covariances. 
for i = 1:101
    ubd = 6*i;
    lbd = ubd - 5;

    A = As(lbd:ubd,1:6);

    %GeqCovProp(:,:.i) = 

end

CeqCovProp = A*CeqCov/A; 



figure; hold on;
ylabels = {'\nu', 'p_1', 'p_2', 'L', 'q_1', 'q_2'}; 

for i = 1:6
subplot(3,2,i); 
plot(nlIndxGeq(:,i),'b');hold on;
p = plot(nlIndxEq(:,i));
plot_latex(p, 'Node', ylabels(i), '', '', {'GEqOE', 'EqOE'});
ylim([0,1.2])
end




