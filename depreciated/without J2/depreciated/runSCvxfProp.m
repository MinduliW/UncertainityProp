function [dynamics,control,eta,fval,egO,dist,TR,nlIndx] =...
    runSCvxfProp(x0Osc, xfMeanKep, xFprop, N, leg1, param)

% FUNCTION NAME:
%   runSCvx
%
% DESCRIPTION:
%   Runs the SCVx algorithm to obtain successive convexification outcomes
%   Convergence based on forward propagation
%
% INPUT:
%   N - (double) number of nodes
%   xFprop - (double []) initial guess of the trajectory
%   param - (struct) Problem parameters
%

if param.coordSys == 1
    RV2Eq = @RV2GEq;
    Eq2RV = @GEq2RV; 
else
    RV2Eq = @RV2CEq;
    Eq2RV = @CEq2RV; 
    
end


slackxf = true;
param.tol = 1e-13;
if slackxf == true
    Scvxfunc = @CoptSolveSlackxf; 
else
    Scvxfunc = @CoptSolvewithoutslack;
end

param.x0OscECI = x0Osc;
param.x0OscEQ(1:6) = RV2Eq(param.x0OscECI(1:3),param.x0OscECI(4:6),param.mu, param.J2, param.Re);

param.xfMeanKep = xfMeanKep;
[rr, vv] = CoordConv.po2pv(xfMeanKep,param.mu);
param.xfMeanCart =  [rr;vv];

% xgO = xFprop;
xguess = xFprop;
param.changenu  =false;
plotres = false;

if plotres == true
    figure;hold on;
    p = plot3(x0Osc(1), x0Osc(2), x0Osc(3), 'r*');
   % p = plot3(xFprop(:,1), xFprop(:,2),xFprop(:,3));

    legnd{1} = 'Departure';
    legnd{2} = 'Target';
    legnd{3} = 'Initial guess';
    plot3(param.xfMeanCart(1), param.xfMeanCart(2), param.xfMeanCart(3), 'b*');

end

% xguesspath = xguess(:,1:3);
% convert xguess to the correct coordinates too
xguessnew  = zeros(size(xguess));
for i = 1: N+1
    xguessnew(i,:) = [RV2Eq(xguess(i,1:3),xguess(i,4:6),param.mu, param.J2, param.Re), xguess(i,7:9)];
end
xguess = xguessnew;

xgO = xguess;



egO = zeros(size(param.tvec));
for val = 1: N+1
    egO(val) =  norm(xFprop(val,7:9));
end

diter = 0;
dst =  leg1.a(1);
param.adddv = 0;
dist = [];
fval = [];
paramArray = createparamArray(param,N);

virtualcxfmax = 100;

dynamicsNew  = zeros(6,N+1);

while diter <10 && (dst > 1e3/param.LU) 

    diter = diter +1;

    legnd{diter+3} = strcat('Iteration:', num2str(diter));
    [As,Bs, constants, TR, xfconv,nlIndx] = runDriver(paramArray, xguess);

    [x,param.changenu] = Scvxfunc(As,Bs,constants,xguess,xgO, N,TR,xfconv, param);

    chanceiters = 0;
    while   param.changenu  == true 

        chanceiters = chanceiters+1;
        param.adddv = param.adddv + 1/N*param.TU/param.LU;

        [As,Bs, constants, TR, xfconv] = runDriver(paramArray, xguess);

        [x,param.changenu] = Scvxfunc(As,Bs,constants,xguess,xgO, N,TR,xfconv, param);
         diter = diter +1;

    end

    dynamics = reshape(x(1:6*(N+1)), [6, N+1]);
    control = reshape( x(6*(N+1)+1: 6*(N+1)+3*N), [3, N]);
    control = [control'; zeros(3,1)'];
    eta = x(6*(N+1)+3*N+1:6*(N+1)+3*N+N);

    fval(diter) = sum(eta);


    for i = 1: N+1
        [R,V] =  Eq2RV(dynamics(:,i),param.tol,param.mu, param.J2, param.Re);
        dynamicsNew(:,i) = [R,V];
    end
    dynamics = dynamicsNew;


if plotres == true
    p = plot3(dynamics(1,1:end-1), dynamics(2,1:end-1),dynamics(3,1:end-1),'+-');
end

    % forward propagate.
    [~,~,~,xguessRV] = fropguess(N, leg1,  x0Osc,param.tvec(end), control,param);

    for i = 1: N+1
        xguessnew(i,:) = [RV2Eq(xguessRV(i,1:3),xguessRV(i,4:6),param.mu, param.J2, param.Re), xguessRV(i,7:9)];
    end
    xguess = xguessnew;


    %% Find difference between current and guess

    % setting dst to be max distance difference between two nodes 
    xnewpath = dynamics(1:3, :)';
    xfproppath = xguessRV(:,1:3);

    etaguess = 0;
    d = [];
    for val = 1:N+1
        etaguess =  etaguess+ norm(xguess(val,7:9));
        d(val) = norm(xfproppath(val,:) - xnewpath(val,:));
    end

    dst = max(d);
    dist = [dist, dst]; 

end


if plotres == true
    plot_latex(p, 'x(m)', 'y(m)','z (m)', '' ,legnd)

    figure;hold on;
    p = plot3(dynamics(1,1:end-1), dynamics(2,1:end-1),dynamics(3,1:end-1),'+-');
    plot3(param.x0OscECI(1), param.x0OscECI(2), param.x0OscECI(3), '*');
  
    plot_latex(p, 'x(m)', 'y(m)','z (m)', 'Final solution',{})

    if param.method == 2
        quiver3(dynamics(1,:),dynamics(2,:),dynamics(3,:),control(:,1)', control(:,2)', control(:,3)',1, 'LineWidth', 2); % plot the delta-Vs
    end
end



end
