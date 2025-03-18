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

ppGeo;
%[rr, vv] = CoordConv.po2pv([a0, e0, inc0, RAAN0_t0,0.0001,0.001], param.mu);
%[rr, vv] = CoordConv.po2pv([a0,  0.0064, deg2rad(99.2248),  deg2rad(20), deg2rad(20), deg2rad(20)],param.mu);

[rr, vv] = CoordConv.po2pv([a0,  0.0064, deg2rad(0.01),  deg2rad(5), deg2rad(5), deg2rad(5)],param.mu);

x0cartOsc = [rr;vv];

nsamples = 100;

% variance 
r_var =  1e3/param.LU ;
v_var =  1/param.LU*param.TU; 

CartX = [];

CartX(1,:) =  x0cartOsc';

mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/propmain3.cpp -ldace -lgsl -lfmt -lcblas

% orbital period initial
periodInit = 2*pi*sqrt(a0^3/param.mu);
param.t0 = 0; 
param.tf = 15*periodInit;
param.nu = 1;
param.method = 2; 

param.coordSys = 1; 
paramArray = createparamArray(param,nsamples+1);


MEEX0   =  CoordConv.vec2mee(CartX(1,1:3),CartX(1,4:6),param.mu);
GEqoeX0 = RV2GEq(CartX(1,1:3),CartX(1,4:6),param.mu, param.J2, param.Re);
CEqoeX0 = RV2CEq(CartX(1,1:3),CartX(1,4:6),param.mu, param.J2, param.Re);
kepx0   = vec2orbElem(CartX(1,1:3),CartX(1,4:6),param.mu);

for i = 2:nsamples+1
    
    val = [randn(1,3),0,0,0]'.*[r_var*ones(3,1); zeros(3,1)] + [0,0,0,randn(1,3)]'.*[zeros(3,1); v_var*ones(3,1)];
    
    dCartX(i,:) =  val'; 
    
    [fxe, fx0 , STM0, STMe, fxMEE, fx0MEE, STM0MEE, STMeMEE,...
        fxCeq, fx0Ceq, STM0Ceq, STMeCeq,  fxGeq, fx0Geq, STM0Geq, STMeGeq, ...
     fxKep(:,i) ,fx0Kep, STM0Kep, STMeKep]= propmain3(paramArray, CartX, dCartX(i,:));
    
    dCartXX(i,:) = fxe - fx0;
    dMEEXX(i,:)  = fxMEE - fx0MEE;
    dCeqXX(i,:)  = fxCeq - fx0Ceq; 
    dGeqXX(i,:)  = fxGeq - fx0Geq; 
    dKepxx(i,:)  = fxKep(:,i) - fx0Kep; 

    nlc(i)    = norm(STMe - STM0)/norm(STM0); 
    nlcMEE(i) = norm(STMeMEE - STM0MEE)/norm(STM0MEE); 
    nlcCeq(i) = norm(STMeCeq - STM0Ceq)/norm(STM0Ceq); 
    nlcGeq(i) = norm(STMeGeq - STM0Geq)/norm(STM0Geq); 
    nlcKep(i) = norm(STMeKep - STM0Kep)/norm(STM0Kep); 


end

%% 

% get the initial error in different coordinates.
CartX0 = CartX(1,:) +dCartX;

for i = 2: nsamples+1

    MEEX1   =  CoordConv.vec2mee(CartX0(i,1:3),CartX0(i,4:6),param.mu);
   
    GEqoeX1 = RV2GEq(CartX0(i,1:3),CartX0(i,4:6),param.mu, param.J2, param.Re);
    CEqoeX1 = RV2CEq(CartX0(i,1:3),CartX0(i,4:6),param.mu, param.J2, param.Re);
    kepx1   = vec2orbElem(CartX0(i,1:3),CartX0(i,4:6),param.mu);

    
    dMEEX0(i,:)   = MEEX1 - MEEX0; 


    dGEqoeX0(i,:) = GEqoeX1 - GEqoeX0; 
    dCEqoeX0(i,:) = CEqoeX1 - CEqoeX0; 
    dkepx0(i,:)   = kepx1 - kepx0; 


end



%% 
figure; 

subplot(1,5,1); hold on; 
 p = plot(dCartXX(:,1), dCartXX(:,2), 'r.');

plot(dCartX(:,1), dCartX(:,2), 'k.');
plot_latex(p, 'X[Re]', 'Y[re]','', 'Cartesian' ,{});


subplot(1,5,2); hold on; 
p = plot(dMEEXX(:,1), dMEEXX(:,6), 'r.');
plot(dMEEX0(:,1), dMEEX0(:,6), 'k.');
plot_latex(p, 'p [Re]', 'L[rad]','', 'MEE' ,{});


subplot(1,5,3); hold on; 
p = plot(dGeqXX(:,1), wrapToPi(dGeqXX(:,4)), 'r.');
plot(dGEqoeX0(:,1), wrapToPi(dGEqoeX0(:,4)), 'k.');
plot_latex(p, '\nu', 'L[rad]','', 'Geqoe' ,{});


subplot(1,5,4); hold on; 
p = plot(dCeqXX(:,1), wrapToPi(dCeqXX(:,4)), 'r.');
plot(dCEqoeX0(:,1), wrapToPi(dCEqoeX0(:,4)), 'k.');
plot_latex(p, '\nu', 'L[rad]','', 'Ceqoe' ,{});



subplot(1,5,5); hold on; 
p = plot(dKepxx(:,1), wrapToPi(dKepxx(:,4)+dKepxx(:,5)+dKepxx(:,6)), 'r.');
plot(dkepx0(:,1), wrapToPi(dkepx0(:,4)+dkepx0(:,5)+dkepx0(:,6)), 'k.');
%plot(dCartX(:,1), dCartX(:,6), 'k.');
plot_latex(p, 'a[Re]', 'L[rad]','', 'Kepler' ,{});


%%

figure; 

subplot(1,5,1); 
semilogy(1:nsamples, nlc(2:end))
hold on; ylim([10e-9,100]);
semilogy(1:nsamples, max(nlc(2:end))*ones(size(nlc(2:end))));
plot_latex(p, 'sample', 'nonlinear index','', 'Cartesian' ,{});

subplot(1,5,2); 
semilogy(1:nsamples,nlcMEE(2:end));
hold on; ylim([10e-9,100]);
semilogy(1:nsamples,max(nlcMEE(2:end))*ones(size(nlcMEE(2:end))));

plot_latex(p, 'sample', 'nonlinear index','', 'MEE' ,{});


subplot(1,5,3); 
semilogy(1:nsamples,nlcGeq(2:end));
hold on; ylim([10e-9,100]);
semilogy(1:nsamples,max(nlcGeq(2:end))*ones(size(nlcGeq(2:end))));
plot_latex(p, 'sample', 'nonlinear index','', 'Geqoe' ,{});

subplot(1,5,4); 
semilogy(1:nsamples,nlcCeq(2:end));
hold on; ylim([10e-9,100]);
semilogy(1:nsamples,max(nlcCeq(2:end))*ones(size(nlcCeq(2:end))));
plot_latex(p, 'sample', 'nonlinear index','', 'Ceqoe' ,{});


subplot(1,5,5); 
semilogy(1:nsamples,nlcKep(2:end));
hold on; 
ylim([10e-9,100]);
semilogy(1:nsamples,max(nlcKep(2:end))*ones(size(nlcKep(2:end))));
plot_latex(p, 'sample', 'nonlinear index','', 'Kepler' ,{});

%% 

figure;
h1 = semilogy(1:nsamples, nlc(2:end),'b')
hold on; 
%ylim([10e-9,10]);
semilogy(1:nsamples, max(nlc(2:end))*ones(size(nlc(2:end))),'b--');

 
h2 = semilogy(1:nsamples,nlcMEE(2:end), 'g');
semilogy(1:nsamples,max(nlcMEE(2:end))*ones(size(nlcMEE(2:end))),'g--');


h3 = semilogy(1:nsamples,nlcGeq(2:end),'r');
semilogy(1:nsamples,max(nlcGeq(2:end))*ones(size(nlcGeq(2:end))),'r--');

 
h4 = semilogy(1:nsamples,nlcCeq(2:end),'c');
semilogy(1:nsamples,max(nlcCeq(2:end))*ones(size(nlcCeq(2:end))),'c--');


h5 = semilogy(1:nsamples,nlcKep(2:end),'m');
semilogy(1:nsamples,max(nlcKep(2:end))*ones(size(nlcKep(2:end))),'m--');
ylim([10^-10, 10^0])
plot_latex(p, 'sample', 'nonlinear index','', '' ,{});

legend([h1 h2, h3, h4, h5],{'Cartesian','MEE', 'GEqoOE', 'CEqoe', 'Keplerian'})

