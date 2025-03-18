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


%problemParam;
ppGeo;
%[rr, vv] = CoordConv.po2pv([a0, e0, inc0, RAAN0_t0,0.0001,0.001], param.mu);
[rr, vv] = CoordConv.po2pv([a0, 0.000001, 0.000001, 0.001,0.0001,0.001],param.mu);
x0cartOsc = [rr;vv];

%x0cartOsc= [ a0, 0, 0, 0.0001, v0, 0.00001]';

nsamples = 100;

% variance 
r_var =  1e2/param.LU ;
v_var =  0.1/param.LU*param.TU; 

CartX = [];

CartX(1,:) =  x0cartOsc';
for i = 2:nsamples+1
    val = x0cartOsc + [randn(1,3),0,0,0]'.*[r_var*ones(3,1); zeros(3,1)] + [0,0,0,randn(1,3)]'.*[zeros(3,1); v_var*ones(3,1)];
    CartX(i,:) =  val'; 
end

% convert each coordinate to Eq.
GEqoeX = [];CEqoeX = [];MEEx0 = [];kepx0 = [];
for i = 1:nsamples+1
    GEqoeX(i,:) = [RV2GEq(CartX(i,1:3),CartX(i,4:6),param.mu, param.J2, param.Re)];
    CEqoeX(i,:) = [RV2CEq(CartX(i,1:3),CartX(i,4:6),param.mu, param.J2, param.Re)];
    MEEx0(i,:) = CoordConv.vec2mee(CartX(i,1:3),CartX(i,4:6),param.mu);
    kepx0(i,:) = vec2orbElem(CartX(i,1:3),CartX(i,4:6),param.mu);

end

%% transfrom to Geqoe. 
mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/propmain.cpp -ldace -lgsl -lfmt -lcblas

% orbital period initial
periodInit = 2*pi*sqrt(a0^3/param.mu);
param.t0 = 0; 
param.tf = 15*periodInit;
param.nu = 1;
param.method = 2; 

param.coordSys = 1; 
paramArray = createparamArray(param,nsamples+1);

[cartxp, cartnl, MEE, MEEnl , Geqoe, Geqoenl, Ceqoe, Ceqoenl, Kep, kepnl ]...
    = propmain(paramArray, CartX);

%%  Uncertainity propagation
figure; 

subplot(1,5,1); hold on; 
for i = 1:nsamples+1
    p = plot(cartxp(:,1)-cartxp(1,1), cartxp(:,2)-cartxp(1,2), 'r.');
end
plot(CartX(:,1)-CartX(1,1), CartX(:,2)-CartX(1,2), 'k.');
plot_latex(p, 'X[Re]', 'Y[re]','', 'Cartesian' ,{});

subplot(1,5,2); hold on;
for i = 1:nsamples+1
    p = plot((MEE(:,1)-MEE(1,1)), ((MEE(:,6)-MEE(1,6))), 'r.');
end
plot((MEEx0(:,1)-MEEx0(1,1)),((MEEx0(:,6)-MEEx0(1,6))), 'k.');
plot_latex(p, 'p [Re]', 'L[rad]','', 'MEE' ,{});


subplot(1,5,3); hold on;
for i = 1:nsamples+1
    p = plot((Geqoe(:,1)-Geqoe(1,1)), wrapToPi((Geqoe(:,4)-Geqoe(1,4))), 'r.');
end
plot((GEqoeX(:,1)-GEqoeX(1,1)),wrapToPi((GEqoeX(:,4)-GEqoeX(1,4))), 'k.');
plot_latex(p, '\nu', 'L[rad]','', 'Geqoe' ,{});


subplot(1,5,4); hold on;
for i = 1:nsamples+1
    p = plot((Ceqoe(:,1)-Ceqoe(1,1)), wrapToPi((Ceqoe(:,4)-Ceqoe(1,4))), 'r.');
end
plot((CEqoeX(:,1)-CEqoeX(1,1)),wrapToPi((CEqoeX(:,4)-CEqoeX(1,4))), 'k.');
plot_latex(p, '\nu', 'L[rad]','', 'Ceqoe' ,{});


subplot(1,5,5); hold on;
truelon = wrapTo2Pi(wrapTo2Pi(Kep(:,6)) +  wrapTo2Pi(Kep(:,5)) +wrapTo2Pi(Kep(:,4))); 

for i = 1:nsamples+1
    p = plot((Kep(:,1)-Kep(1,1)), wrapToPi((truelon-truelon(1))), 'r.');
end
plot((kepx0(:,1)-kepx0(1,1)),wrapToPi((kepx0(:,7)-kepx0(1,7))), 'k.');
plot_latex(p, 'a [Re]', 'L [rad]','', 'Kepler' ,{});

% 
% figure; hold on; 
% plot((wrapToPi(Kep(:,6)-Kep(1,6)))); 
% plot((wrapToPi(kepx0(:,6)-kepx0(1,6))),'--'); 
% 
% 

%% Nonlinearity index. 

figure; 

subplot(1,5,1); hold on; 
for i = 1
    plot(cartnl(:,i));
end
plot_latex(p, 'sample', 'nonlinear index','', 'Cartesian' ,{'x', 'y', 'z', 'v_x', 'v_y', 'v_z'});


subplot(1,5,2); hold on; 
for i = 1
    plot(MEEnl(:,i));
end
plot_latex(p, 'sample', 'nonlinear index','', 'MEE' ,{'p', 'f', 'g', 'h', 'k', 'L'});



subplot(1,5,3); hold on; 
for i = 1
    plot(Geqoenl(:,i));
end
plot_latex(p, 'sample', 'nonlinear index','', 'Geqoe' ,{'\nu', 'p_1', 'p_w', 'L', 'q_1', 'q_2'});



subplot(1,5,4); hold on; 
for i = 1
    plot(Ceqoenl(:,i));
end
plot_latex(p, 'sample', 'nonlinear index','', 'Ceqoe' ,{'\nu', 'p_1', 'p_w', 'L', 'q_1', 'q_2'});


subplot(1,5,5); hold on; 
for i = 1
    plot(kepnl(:,i));
end
plot_latex(p, 'sample', 'nonlinear index','', 'Kepler' ,{'a', 'e', 'i', '\Omega', '\omega', '\theta'});



max(cartnl)
max(MEEnl)
max(Geqoenl)
max(Ceqoenl)
max(kepnl)