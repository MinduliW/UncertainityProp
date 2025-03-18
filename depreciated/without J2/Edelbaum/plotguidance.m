function plotguidance(x0, a0, af, inc0, incf, RAAN0, lastOmega, period, betaedel,f,t2,param)


param.eccAdj = true ;

etas = 1; %0:0.2:1;


param.etae = 0;
param.etaa = 0;
param.etai = 0;


as = [];
es = [];
incs = [];
RAANs = [];

for i = 1: length(etas)
    param.eps = etas(i);
    
   
    [time_ga,states_ga] = ode45(@(t,x) ...
        forwardPropagator(t,x, betaedel,f, t2, param),...
        [0:period/20:t2(end)],[x0';param.m0]);
    fprintf('Fuel burnt Edelbaum + Forward Propagation +Guidance adjustments : %f kg \n',...
    param.m0 - states_ga(end,end));


    dv = param.Isp*param.g0*log(param.m0/states_ga(end,end))
    for j = 1 : length(time_ga)
    [~, ~, ~, ~, ~, coemean(j,:),~] =...
        forwardPropagator(time_ga(j),states_ga(j,1:7)',...
        betaedel,f , t2, param);
    end
    
    time(i,:) = time_ga;
    as(i,:) = coemean(:,1);
    es(i,:) = coemean(:,2);
    incs(i,:) =  coemean(:,3); 
    RAANs(i,:) = coemean(:,4); 
end


figure; 
subplot(1,4,1);hold on;

plot(time_ga/param.TU, a0*ones(size(time_ga))/1e3, 'k--');
plot(time_ga/param.TU, af*ones(size(time_ga))/1e3, 'g--');
for i = 1: length(etas)
p = plot(time(i,:)/param.TU, as(i,:)/1e3);
end
plot_latex(p ,'time (days)', 'sma (km)', '', '', {'Initial', 'Final', ...
    '$\epsilon = 0$','$\epsilon = 0.2$','$\epsilon = 0.4$','$\epsilon = 0.6$',...
    '$\epsilon = 0.8$', '$\epsilon = 1$' })

subplot(1,4,2); hold on;
plot(time_ga/param.TU, inc0*ones(size(time_ga))*180/pi, 'k--');
plot(time_ga/param.TU, incf*ones(size(time_ga))*180/pi, 'g--');
for i = 1: length(etas)
p = plot(time(i,:)/param.TU, incs(i,:)*180/pi);
end
plot_latex(p ,'time (days)', 'inc (deg)', '', '', {})


subplot(1,4,3);hold on;
plot(t2/param.TU, RAAN0*ones(size(t2))*180/pi, 'k--');
plot(t2/param.TU, lastOmega*ones(size(t2))*180/pi, 'g--');
for i = 1: length(etas)
p = plot(time(i,:)/param.TU, unwrap(RAANs(i,:))*180/pi);
end

plot_latex(p ,'time (days)', 'Omega (deg)', '', '', {})

subplot(1,4,4); hold on; 
for i = 1: length(etas)
p = plot(time(i,:)/param.TU, es(i,:));
end

plot_latex(p ,'time (days)', 'eccentricity', '', '', {})



    
end


