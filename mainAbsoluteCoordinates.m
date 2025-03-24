clear,clc; close all;

% This is a cleaned up version of the nonlinearity index calculation.
% What is being done:
    % 

%% Initial orbit 
mu = 3.986005e14;%m^2/s^3
Re = 6378.137e3; % m
J2 =  1.08262668e-3;
g0 = 9.80665; %m/s
m_spacecraft = 800;

% normalization of units for easier integration. 
param.LU = Re;
param.VU = sqrt(mu/Re);
param.TU = param.LU/param.VU;
param.MU = m_spacecraft;

param.mu= mu/param.LU^3*param.TU^2; 
param.Re = Re/param.LU;
param.J2 = J2; 

a0 = 350e3/param.LU + param.Re;

% convert initial state to cartesian. 
[rr, vv] = CoordConv.po2pv([a0,  0.0064, deg2rad(0.01),  deg2rad(5), ...
    deg2rad(5), deg2rad(5)],param.mu);

CartX(1,:) = [rr;vv];

% orbital period initial
periodInit = 2*pi*sqrt(a0^3/param.mu);
param.t0 = 0; 
param.tf = 15*periodInit;

paramArray = [param.mu, param.J2, param.Re, param.tf];

% Get the initial state in different coordinates 
MEEX0 =  CoordConv.vec2mee(CartX(1:3),CartX(4:6),param.mu);
GEqoeX0 = CoordConv.RV2GEq(CartX(1:3),CartX(4:6),param.mu, param.J2,param.Re);
CEqoeX0 = CoordConv.RV2CEq(CartX(1:3),CartX(4:6),param.mu, param.J2,param.Re);
kepx0   = CoordConv.vec2orbElem(CartX(1:3),CartX(4:6),param.mu);

%% Propagate nominal trajectory 
mex -v -R2018a -I/usr/local/include -L/usr/local/lib cpp/CloudProp.cpp -ldace 

% this function propagates the spacecraft by 15 orbits, and also gives you
% the STMs. fx0s are the final states with different coordinate systems,
% and STM0s are the STMs in different coordinate systems 
[fx02, STM0,fx0MEE, STM0MEE, fx0Geq, STM0Geq, fx0Ceq, STM0Ceq, fx0kep, STM0kep] = ...
    CloudProp(paramArray, CartX,MEEX0, GEqoeX0, CEqoeX0, kepx0);


%% Perturbed state propagation 

nPoints = 50; % number of pertubed states we are looking at. 

% The covariance of the pertubation introduced in cartesian coordinates. It
% is 1 km and 1 m/s here. 
r_var = 1e3/param.LU ;
v_var = 1/param.LU*param.TU; 


for i = 1: nPoints

    val(i,:) = [randn(1,3), 0,0,0]'.*[r_var*ones(3,1); zeros(3,1)] + ...
        [0,0,0,randn(1,3)]'.*[zeros(3,1); v_var*ones(3,1)];

    % Update the cartesian initial state to the perturbed initial state 
    CartXnew = CartX+val(i,:); 

    % Convertto different coordinate systems 
    MEEX =  CoordConv.vec2mee(CartXnew(1:3),CartXnew(4:6),param.mu);
    GEqoeX = CoordConv.RV2GEq(CartXnew(1:3),CartXnew(4:6),param.mu, param.J2,param.Re);
    CEqoeX = CoordConv.RV2CEq(CartXnew(1:3),CartXnew(4:6),param.mu, param.J2,param.Re);
    kepX  = CoordConv.vec2orbElem(CartXnew(1:3),CartXnew(4:6),param.mu);

    % In each coordinate system, whats the difference between the perturbed
    % initial state and nominal initial state? 
    valMEE(i,:)  = MEEX- MEEX0; 
    valGeqOE(i,:) = GEqoeX- GEqoeX0; 
    valCeqOE(i,:) = CEqoeX- CEqoeX0; 
    valkep(i,:)   = kepX - kepx0;

    % Like before, we propagate the now perturbed initial state to get the
    % perturbed state transition matricies and the perturbed final states.
    [fxe2, STMe,fxeMEE, STMeMEE, fxeGeq, STMeGeq,fxeCeq, STMeCeq , fxekep(i,:), STMekep] =...
        CloudProp(paramArray, CartXnew,MEEX,GEqoeX, CEqoeX, kepX);

    % This is a calculaton of the nonlineartiy index based on the frobenius
    % norm from Eq. 5 in
    % https://link.springer.com/article/10.1007/BF03546420.
    nuCart(i) = norm(STMe-STM0,'fro')/norm(STM0,'fro');
    nuMEE(i) = norm(STMeMEE-STM0MEE,'fro')/norm(STM0MEE,'fro');
    nuGeq(i) = norm(STMeGeq-STM0Geq,'fro')/norm(STM0Geq,'fro');
    nuCeq(i) = norm(STMeCeq-STM0Ceq,'fro')/norm(STM0Ceq,'fro');
    nuKep(i) = norm(STMekep-STM0kep,'fro')/norm(STM0kep,'fro');

    % Now we calculate how much the perturbed state differes from the
    % nominal state after propagation in each coordinate system. 
    dCartXXnlc(i,:) = fxe2 - fx02; 
    dMEEXXnlc(i,:) = fxeMEE - fx0MEE;
    dGeqXXnlc(i,:) = fxeGeq - fx0Geq;
    dCeqXXnlc(i,:) = fxeCeq - fx0Ceq;
    dkepXXnlc(i,:) = fxekep(i,:)' - fx0kep;
    
end

%% Analysis 

% This plots the nonlinearity indicies - this is what you must compare for
% relative coordinate systems  
figure;  
loglog(max(nuCart), 'o');hold on; 
loglog(max(nuMEE), 'o');
loglog(max(nuGeq), 'o');
loglog(max(nuGeq), 'o');
p = loglog(max(nuKep), 'o');
ylim([1e-7,1e2])
plot_latex(p, '', '\eta','', 'Nonlinearity index from point cloud analysis', ...
    {'Cartesian', 'MEE', 'Geqoe', 'Ceqoe', 'Keplerian'});

% This is a distribution of uncertainity - a linear trend here is good,
% obviously the flatter the red curve is the better/less uncertainity as
% well. 
figure; 
subplot(1,5,1);hold on; 
p = plot(dCartXXnlc(:,1), dCartXXnlc(:,2), 'r.');
% p = plot(dCartXX(:,1), dCartXX(:,2), 'b.');
plot(val(:,1), val(:,2), 'k.');
plot_latex(p, 'X[Re]', 'Y[re]','', 'Cartesian' ,{'Cloud prop', 'Initial' });

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
