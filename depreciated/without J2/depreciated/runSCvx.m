function [dynamics,control,eta,fval,egO] =...
    runSCvx(x0Osc, xfMeanKep, xguess, N, param)

% FUNCTION NAME:
%   runSCvx
%
% DESCRIPTION:
%   Runs the SCVx algorithm to obtain successive convexification outcomes
%
% INPUT:
%   N - (double) number of nodes 
%   xguess - (double []) initial guess of the trajectory
%   param - (struct) Problem parameters 
%

param.x0OscECI = x0Osc;
param.xfMeanKep = xfMeanKep;
[rr, vv] = CoordConv.po2pv(xfMeanKep,param.mu);
param.xfMeanCart =  [rr;vv];

param.changenu =false;

param.x0OscEQ(1:6) = RV2GEq(param.x0OscECI(1:3),param.x0OscECI(4:6),param.mu,param.J2, param.Re);

xguesspath = xguess(:,1:3);
% convert xguess to the correct coordinates too
for i = 1: length(xguess)
    xguessnew(i,:) = [RV2GEq(xguess(i,1:3),xguess(i,4:6),param.mu,param.J2, param.Re), xguess(i,7:9)];
end
xguess = xguessnew;




xgO = xguess;

fileID = fopen('iter.txt','w');
fprintf(fileID,'%12.49f %12.49f %12.49f %12.49f %12.49f %12.49f %12.49f %12.49f %12.49f\n',xguess');
fclose(fileID);

plotres = false;

if plotres == true
    figure;hold on;
    plot3(x0Osc(1), x0Osc(2), x0Osc(3), 'r*');
    p = plot3(xguess(:,1), xguess(:,2),xguess(:,3));

    legnd{1} = 'Departure';
    legnd{2} = 'Target';
    legnd{3} = 'Initial guess';
    plot3(param.xfMeanCart(1), param.xfMeanCart(2), param.xfMeanCart(3), 'b*');

end

egO = zeros(size(param.tvec));
for val = 1: length(xguess)
    egO(val) =  norm(xguess(val,7:9));
end

    

distance = 100;diter = 0;

param.adddv = 0; 

if ismac
    ! g++ -o main main.cpp -I/usr/local/include -L/usr/local/lib -ldace -lgsl -std=c++11

else
    !wsl g++ -o main main.cpp -I/usr/local/include -L/usr/local/lib -ldace -lgsl -std=c++11

end

while diter <20  % && distance >0.05
    if diter > 0
        xguesspath = xnewpath;
    end
    diter = diter +1;
    legnd{diter+3} = strcat('Iteration:', num2str(diter));

    if ismac
        !./main
    else
        ! wsl ./main
    end

    %%

    [As, Bs,constants,TR,xfconv] =extractfromtext('outputcpp.txt','trustRegion.txt', 'xfconv.txt');

    [x,param.changenu] = CoptSolvewithoutslack(As,Bs,constants,xguess,xgO, N,TR,xfconv, param);


    chanceiters = 0;
    while   param.changenu  == true % in case the trust region is too small

        chanceiters = chanceiters+1;
        param.adddv = param.adddv + 1/N*param.TU/param.LU;

        fileID = fopen('iter.txt','w');
        fprintf(fileID,'%12.49f %12.49f %12.49f %12.49f %12.49f %12.49f %12.49f %12.49f %12.49f\n',xgO');
        fclose(fileID);
        if ismac
            !./runmac
        else
            !wsl ./main
        end

        [As, Bs,constants,TR,xfconv] =extractfromtext('outputcpp.txt','trustRegion.txt', 'xfconv.txt');
        [x,param.changenu] = CoptSolvewithoutslack(As,Bs,constants,xgO,xgO, N,TR,xfconv,param);

    end

    dynamics = reshape(x(1:6*(N+1)), [6, N+1]);
    control = reshape( x(6*(N+1)+1: 6*(N+1)+3*N), [3, N]);
    control = [control'; zeros(3,1)'];
    eta = x(6*(N+1)+3*N+1:6*(N+1)+3*N+N);

    fval(diter) = sum(eta);

    xnew = [dynamics',control];

    for i = 1: length(xguess)
        [R,V] = GEq2RV(dynamics(:,i),1e-13, param.mu,param.J2, param.Re);
        dynamicsNew(i,:) = [R,V];
    end
    dynamics = dynamicsNew';

    if plotres == true
        p = plot3(dynamics(1,1:end-1), dynamics(2,1:end-1),dynamics(3,1:end-1), '+-');
    end


    % write xguess into a text file.
    fileID = fopen('iter.txt','w');
    fprintf(fileID,'%12.49f %12.49f %12.49f %12.49f %12.49f %12.49f %12.49f %12.49f %12.49f\n',xnew');
    fclose(fileID);

    %% Find difference between current and guess

    % xguesspath = xguess(:,1:3);
    xnewpath = dynamics(1:3,:)';

    etaguess = 0;
    d = [];
    for val = 1:length(xnewpath)-1
        etaguess =  etaguess+ norm(xguess(val,7:9));
        d(val) = norm(xguesspath(val,:) - xnewpath(val,:));
    end
    distance = sum(d);

    if plotres == true
        figure; subplot(2,1,1);
        plot((xguess(:,7:9)));

        subplot(2,1,2);
        plot((xnew(:,7:9)));
    end

    xguess = xnew;


    fprintf("total dv = %f m/s \n", sum(eta)*param.LU/param.TU) ;

    fprintf("total fprop prev= %f m/s \n", sum(etaguess)*param.LU/param.TU) ;


    distance = abs(sum(etaguess)*param.LU/param.TU - sum(eta)*param.LU/param.TU);

    fprintf("total dv difference = %f m/s \n", distance) ;




end


if plotres == true
    plot_latex(p, 'x(m)', 'y(m)','z (m)', '' ,legnd)

    figure;hold on;
    p = plot3(dynamics(1,1:end-1), dynamics(2,1:end-1),dynamics(3,1:end-1),'+-');
    plot3(param.x0OscECI(1), param.x0OscECI(2), param.x0OscECI(3), '*');
    %plot3(param.xfMeanKep(1), param.xfMeanKep(2), param.xfMeanKep(3), '*');

    plot_latex(p, 'x(m)', 'y(m)','z (m)', 'Final solution',{})

    if param.method == 2
        quiver3(dynamics(1,:),dynamics(2,:),dynamics(3,:),control(:,1)', control(:,2)', control(:,3)',1, 'LineWidth', 2); % plot the delta-Vs
    end
end



end
