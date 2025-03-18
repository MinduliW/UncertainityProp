function [tProp,x0cartOsc,coemean,xguess2] = frop(N, leg1,  x0cartOsc,tendRun, control,param)

% apply the control we get to param.nOrbitsRun, or till tendRun
% apply the control to only

rs = x0cartOsc;

% create an array till tendrun.
tProp = linspace(param.t0,tendRun, N+1);

statesAfterdv = zeros(N+1,6);
coemean =  zeros(N+1,6);

coe = CoordConv.vec2orbElem(rs(1:3),rs(4:6),param.mu);

coemean(1,:) = osc2meankep(coe, param.mu, param.J2, param.Re);

xguess2 =[];

for i = 1: N

    controlval(1) = interp1(param.tvec, control(:,1), tProp(i));
    controlval(2) = interp1(param.tvec, control(:,2), tProp(i));
    controlval(3) = interp1(param.tvec, control(:,3), tProp(i));

    res = [rs(1:3)',rs(4:6)',controlval];

    xguess2 = [xguess2; res];

    % add dv to the current state
    rs(4:6) = rs(4:6)+controlval';

    % store the state with added dv
    statesAfterdv(i,:) = rs;

    % propagate from current node to the next.
    [~,propstate] = ode45(@(t,x) propagateCart(t, x, param),...
        [tProp(i),tProp(i+1)],rs, odeset('RelTol', 1e-13,'AbsTol',1e-13));


    % make rs the propagated state.
    rs =propstate(end,:)';

    coe = CoordConv.vec2orbElem(propstate(end,1:3),propstate(end,4:6),param.mu);

   % coe= [a/param.LU, ecc, incl*pi/180, RAAN*pi/180, argp*pi/180, nu*pi/180];
    coemean(i+1,:)  = osc2meankep(coe, param.mu, param.J2, param.Re);

end
x0cartOsc = rs;

resOsc = [rs',0,0,0] ;
xguess2 = [xguess2; resOsc];


% what is the reference at this point ? 
at =  interp1(leg1.t,leg1.a,tendRun);
inct =  interp1(leg1.t,leg1.inc,tendRun);
RAANt = interp1(leg1.t,leg1.RAAN,tendRun);
% 
% fprintf('a error (km) = %f \n', abs(at- coemean(end,1))*param.LU/1e3);
% fprintf('inc error (deg) = %f \n', abs(inct- coemean(end,3))*180/pi);
% fprintf('RAAN error (deg) = %f \n', acos(cos(RAANt- coemean(end,4)))*180/pi);

end