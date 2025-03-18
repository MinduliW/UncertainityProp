%% nonlinear index comparison
clear,clc; close all;

% Note: this is phasing 
% TO DO: solve a standard convex problem using this
addpath('Edelbaum');
addpath('Dynamics');
addpath('Utility');
addpath('SGP4routines_NAIF');

addmosek;
cspice_furnsh('SGP4routines_NAIF/kernelmac.txt')

ppLEO;
%ppGeo; 

[rr, vv] = CoordConv.po2pv([a0,  0.0064, deg2rad(99.2248),  deg2rad(5), ...
    deg2rad(5), deg2rad(5)],param.mu);

% variance 
r_var = 1e3/param.LU ;
v_var = 1/param.LU*param.TU; 

CartX(1,:) = [rr;vv];

% orbital period initial
periodInit = 2*pi*sqrt(a0^3/param.mu);
param.t0 = 0; 
param.tf = 15*periodInit;
paramArray = createparamArray(param);

dxCartMax = 3*[r_var, r_var,r_var, v_var , v_var, v_var];

MEEX0 =  CoordConv.vec2mee(CartX(1:3),CartX(4:6),param.mu);
GEqoeX0 = RV2GEq(CartX(1:3),CartX(4:6),param.mu, param.J2,param.Re);
CEqoeX0 = RV2CEq(CartX(1:3),CartX(4:6),param.mu, param.J2,param.Re);
kepx0   = vec2orbElem(CartX(1:3),CartX(4:6),param.mu);

mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/CloudProp.cpp -ldace 

[fx02, STM0,fx0MEE, STM0MEE, fx0Geq, STM0Geq, fx0Ceq, STM0Ceq, fx0kep, STM0kep] = ...
    CloudProp(paramArray, CartX,MEEX0, GEqoeX0, CEqoeX0, kepx0);

for i = 1: 100

    val(i,:) = [randn(1,3), 0,0,0]'.*[r_var*ones(3,1); zeros(3,1)] + ...
        [0,0,0,randn(1,3)]'.*[zeros(3,1); v_var*ones(3,1)];

    CartXnew = CartX+val(i,:); 
    MEEX =  CoordConv.vec2mee(CartXnew(1:3),CartXnew(4:6),param.mu);
    GEqoeX = RV2GEq(CartXnew(1:3),CartXnew(4:6),param.mu, param.J2,param.Re);
    CEqoeX = RV2CEq(CartXnew(1:3),CartXnew(4:6),param.mu, param.J2,param.Re);
    kepX  = vec2orbElem(CartXnew(1:3),CartXnew(4:6),param.mu);

    valMEE(i,:)  = MEEX- MEEX0; 
    valGeqOE(i,:) = GEqoeX- GEqoeX0; 
    valCeqOE(i,:) = CEqoeX- CEqoeX0; 
    valkep(i,:)   = kepX - kepx0;

    [fxe2, STMe,fxeMEE, STMeMEE, fxeGeq, STMeGeq,fxeCeq, STMeCeq , fxekep(i,:), STMekep] =...
        CloudProp(paramArray, CartXnew,MEEX,GEqoeX, CEqoeX, kepX);

    dCartXXnlc(i,:) = fxe2 - fx02; 
    dMEEXXnlc(i,:) = fxeMEE - fx0MEE;
    dGeqXXnlc(i,:) = fxeGeq - fx0Geq;
    dCeqXXnlc(i,:) = fxeCeq - fx0Ceq;
    dkepXXnlc(i,:) = fxekep(i,:)' - fx0kep;
    
end

%%

figure; 
subplot(1,5,1);hold on; 
p = plot(dCartXXnlc(:,1), dCartXXnlc(:,2), 'r.');
plot(val(:,1), val(:,2), 'k.');
plot_latex(p, 'X[Re]', 'Y[re]','', 'Cartesian' ,{'Cloud prop','Initial' });

subplot(1,5,2);hold on; 
p = plot(dMEEXXnlc(:,1), dMEEXXnlc(:,6), 'r.');
plot(valMEE(:,1), valMEE(:,6), 'k.');
plot_latex(p, 'p [Re]', 'L[rad]','', 'Modified Equinoctials' ,{});

subplot(1,5,3);hold on; 
p = plot(dGeqXXnlc(:,1), dGeqXXnlc(:,4), 'r.');
plot(valGeqOE(:,1), valGeqOE(:,4), 'k.');
plot_latex(p, '\nu', 'L[rad]','', 'GEqoOE' ,{});

subplot(1,5,4);hold on; 
p = plot(dCeqXXnlc(:,1), dCeqXXnlc(:,4), 'r.');
plot(valCeqOE(:,1), valCeqOE(:,4), 'k.');
plot_latex(p, '\nu', 'L[rad]','', 'CEqoOE' ,{});


subplot(1,5,5); hold on; 
p = plot(dkepXXnlc(:,1), wrapToPi(dkepXXnlc(:,4)+dkepXXnlc(:,5)+dkepXXnlc(:,6)), 'r.');
plot(valkep(:,1),  wrapToPi(valkep(:,4)+valkep(:,5)+valkep(:,6)), 'k.');
plot_latex(p, 'a[Re]', 'L[rad]','', 'Kepler' ,{});

 

