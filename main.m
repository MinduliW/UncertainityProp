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
param.tf = 15*periodInit;
paramArray = createparamArray(param);

dxCartMax = 3*[r_var, r_var,r_var, v_var , v_var, v_var];

MEEX0 =  CoordConv.vec2mee(CartX(1:3),CartX(4:6),param.mu);
GEqoeX0 = RV2GEq(CartX(1:3),CartX(4:6),param.mu, param.J2,param.Re);
CEqoeX0 = RV2CEq(CartX(1:3),CartX(4:6),param.mu, param.J2,param.Re);
kepx0   = vec2orbElem(CartX(1:3),CartX(4:6),param.mu);

% propagate for 15 orbits and get the absol. pos. 
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
[t,y] = ode45(@(t,y) propagateCart(t,y, param), [0, param.tf], CartX, options);

CartXend = y(end,:);

MEEXf   = CoordConv.vec2mee(CartXend(1:3),CartXend(4:6),param.mu);
GEqoeXf = RV2GEq(CartXend(1:3),CartXend(4:6),param.mu, param.J2,param.Re);
CEqoeXf = RV2CEq(CartXend(1:3),CartXend(4:6),param.mu, param.J2,param.Re);
kepxf   = vec2orbElem(CartXend(1:3),CartXend(4:6),param.mu);


% mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/CloudProp.cpp -ldace 
% 
[fx02, STM0,fx0MEE, STM0MEE, fx0Geq, STM0Geq, fx0Ceq, STM0Ceq, fx0kep, STM0kep] = ...
    CloudProp(paramArray, CartX,MEEX0, GEqoeX0, CEqoeX0, kepx0);



%% Covariance method 
% 
% mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/CovarianceProp.cpp -ldace
% 
% for i = 1: 100
% 
%     %valkep(i,:) = kepX - kepx0;
%     val(i,:)  =  1/3*randn(1,6);
% 
%     [fxe, fx0 , nlc,fxeMEE, fx0MEEc, nlcMEE, fxeGeq, fx0Geqc, nlcGeq,...
%        fxeCeq, fx0Ceqc, nlcCeq, fxeKep(i,:), fx0kepc, nlcKep]= ...
%         CovarianceProp(paramArray, CartX , val(i,:));
% 
%     dCartXX(i,:) = fxe- fx0;
%     dMEEXX(i,:) = fxeMEE -fx0MEEc;
%     dGEqoeXX(i,:) = fxeGeq - fx0Geqc;
%     dCEqoeXX(i,:) = fxeCeq - fx0Ceqc;
%     dKepXX(i,:) = fxeKep(i,:)' - fx0kepc;
% end
% 

%% Point cloud method 

for i = 1: 10

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

    nuCart(i) = norm(STMe-STM0,'fro')/norm(STM0,'fro');
    nuMEE(i) = norm(STMeMEE-STM0MEE,'fro')/norm(STM0MEE,'fro');
    nuGeq(i) = norm(STMeGeq-STM0Geq,'fro')/norm(STM0Geq,'fro');
    nuCeq(i) = norm(STMeCeq-STM0Ceq,'fro')/norm(STM0Ceq,'fro');
    nuKep(i) = norm(STMekep-STM0kep,'fro')/norm(STM0kep,'fro');

    dCartXXnlc(i,:) = fxe2 - fx02; 
    dMEEXXnlc(i,:) = fxeMEE - fx0MEE;
    dGeqXXnlc(i,:) = fxeGeq - fx0Geq;
    dCeqXXnlc(i,:) = fxeCeq - fx0Ceq;
    dkepXXnlc(i,:) = fxekep(i,:)' - fx0kep;
    
end

%% 
figure;
% subplot(1,2,1);
% loglog(nlc, 'o'); hold on; 
% loglog(nlcMEE, 'o');
% loglog(nlcGeq, 'o');
% loglog(nlcCeq, 'o');
% p = loglog(nlcKep, 'o');
% ylim([1e-7,1e2])
% plot_latex(p, '', '\eta','', 'Nonlinearity index from Covariance' , ...
%     {'Cartesian', 'MEE', 'Geqoe', 'Ceqoe', 'Keplerian'});

subplot(1,2,2); 
loglog(max(nuCart), 'o');hold on; 
loglog(max(nuMEE), 'o');
loglog(max(nuGeq), 'o');
loglog(max(nuGeq), 'o');
p = loglog(max(nuKep), 'o');
ylim([1e-7,1e2])
plot_latex(p, '', '\eta','', 'Nonlinearity index from point cloud analysis', ...
    {'Cartesian', 'MEE', 'Geqoe', 'Ceqoe', 'Keplerian'});

%% 
figure; 
subplot(1,5,1);hold on; 
p = plot(dCartXXnlc(:,1), dCartXXnlc(:,2), 'r.');
% p = plot(dCartXX(:,1), dCartXX(:,2), 'b.');
plot(val(:,1), val(:,2), 'k.');
plot_latex(p, 'X[Re]', 'Y[re]','', 'Cartesian' ,{'Cloud prop','Cov. Prop', 'Initial' });

subplot(1,5,2);hold on; 
p = plot(dMEEXXnlc(:,1), dMEEXXnlc(:,6), 'r.');
% p = plot(dMEEXX(:,1), wrapToPi(dMEEXX(:,6)), 'b.');
plot(valMEE(:,1), valMEE(:,6), 'k.');
plot_latex(p, 'p [Re]', 'L[rad]','', 'Modified Equinoctials' ,{});

subplot(1,5,3);hold on; 
p = plot(dGeqXXnlc(:,1), dGeqXXnlc(:,4), 'r.');
% p = plot(dGEqoeXX(:,1), wrapToPi(dGEqoeXX(:,4)), 'b.');
plot(valGeqOE(:,1), valGeqOE(:,4), 'k.');
plot_latex(p, '\nu', 'L[rad]','', 'GEqoOE' ,{});

subplot(1,5,4);hold on; 
p = plot(dCeqXXnlc(:,1), dCeqXXnlc(:,4), 'r.');
% p = plot(dCEqoeXX(:,1), wrapToPi(dCEqoeXX(:,4)), 'b.');
plot(valCeqOE(:,1), valCeqOE(:,4), 'k.');
plot_latex(p, '\nu', 'L[rad]','', 'CEqoOE' ,{});


subplot(1,5,5); hold on; 
p = plot(dkepXXnlc(:,1), wrapToPi(dkepXXnlc(:,4)+dkepXXnlc(:,5)+dkepXXnlc(:,6)), 'r.');
% p = plot(dKepXX(:,1), wrapToPi(dKepXX(:,4)+dKepXX(:,5)+dKepXX(:,6)), 'b.');
plot(valkep(:,1),  wrapToPi(valkep(:,4)+valkep(:,5)+valkep(:,6)), 'k.');
plot_latex(p, 'a[Re]', 'L[rad]','', 'Kepler' ,{});


%%

% figure; hold on;
% p = plot(fxekep(:,1), wrapToPi(fxekep(:,6)), 'b.');
% 
% p = plot(fxeKep(:,1), wrapToPi(fxeKep(:,6)), 'r.');
