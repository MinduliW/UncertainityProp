function [x,changenu] = CoptSolvewithoutslack(As,Bs,constants,xg,xgO, N,TR, xfconv, param)

% FUNCTION NAME:
%   CoptSolve
%
% DESCRIPTION:
%   Function runs MOSEK on the given problem
%
% INPUT:
%   As - (double []) STM
%   Bs - (double []) Gamma
%   constants - (double []) Constant part of the state
%   xg - (double []) guess for the trajectory and control
%   N - (double) number of nodes 
%   param - (struct) problem parameters
%
% OUTPUT:
%   x - (double []) Solution (trajectory+control)

clear prob;

paramMosek=[];
paramMosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-14;
paramMosek.MSK_DPAR_INTPNT_TOL_PFEAS = 1e-14;

prob.c   = [zeros(1,6*(N+1)+3*N), ones(1, N)]; %+...

% Dynamics and BCs
Aeq3 = [];Beq3 = [];
for Nnodes = 1: N

    ubd = 6*Nnodes;
    lbd = ubd - 5;

    A = As(lbd:ubd,:);
    B = Bs(lbd:ubd, :);

    dynamicsArray = [zeros(6, 6*(Nnodes-1)), A, -eye(6), zeros(6,6*(N-Nnodes)),...
        zeros(6, 3*(Nnodes-1)), B, zeros(6, 3*(N-Nnodes)), zeros(6,N)];

    constantsarray =  -constants(lbd:ubd) + A*xg(Nnodes,1:6)' + B*xg(Nnodes,7:9)';

    Aeq3 = [Aeq3; dynamicsArray];
    Beq3 = [Beq3; constantsarray];
end

% initial position.
Aeq1 = ([eye(6), zeros(6,10*N)]);
Beq1 = param.x0OscEQ';


Aeq21 =  [zeros(1,6*N), xfconv.A(1,:), zeros(1,4*N)]; 
Beq21 = param.xfMeanKep(1) -xfconv.B(1) +  xfconv.A(1,:)*xg(end,1:6)';

Aeq22 = [zeros(1,6*N), xfconv.A(3,:), zeros(1,4*N)]; 
Beq22 =  param.xfMeanKep(3) -xfconv.B(3)+  xfconv.A(3,:)*xg(end,1:6)';

Aeq23 = [zeros(1,6*N), xfconv.A(4,:), zeros(1,4*N)];
Beq23 = param.xfMeanKep(4) -xfconv.B(4) +  xfconv.A(4,:)*xg(end,1:6)';

Aeq2 = [Aeq21;Aeq22; Aeq23];
Beq2 = [Beq21;Beq22; Beq23];
% 
% Aeq2 =  [zeros(6,6*N), xfconv.A, zeros(6,4*N)]; 
% Beq2 = param.xfMeanKep' - xfconv.B + xfconv.A*xg(end,1:6)';


% 
% Aeq2 = [zeros(6,6*N), eye(6), zeros(6,4*N)];
% Beq2 = param.xfCart;

Aeq = [Aeq1;Aeq2;Aeq3];Beq = [Beq1;Beq2;Beq3];
%Aeq = [Aeq1;Aeq3];Beq = [Beq1;Beq3];


prob.a = sparse(Aeq); % store matrix efficiently as a sparse matrix

% Lower and Upper Equality Constraints i.e. Ax=b
prob.blc = Beq;
prob.buc = Beq;


%compute spacecraft mass

%%
mass(1) = param.m0;
for i = 2:N+1
    
    dv = norm(xg(i-1,7:9));
    mass(i) = mass(i-1)/exp(dv/param.Isp/param.g0);   
end

%% Setting blc and bux by the nonlinear index.
prob.blx  =[];
prob.bux =[];

TR = TR';

for i = 1: N
    prob.bux  = [prob.bux , xg(i,1:3)+TR(1:3,i)',xg(i,4:6)+TR(4:6,i)'];
    prob.blx = [prob.blx, xg(i,1:3)-TR(1:3,i)',xg(i,4:6)-TR(4:6,i)'];
end

% prob.bux  = [prob.bux , 10*ones(6*N,1)'];
% prob.blx = [prob.blx, -10*ones(6*N,1)'];

prob.bux  = [prob.bux , 10*ones(6,1)'];
prob.blx = [prob.blx, -10*ones(6,1)'];



%% control limits 

if param.impulse == true
    m1 = 3;

else
    m1 = 1;

end

% introduce eclipses
if param.eclipses == true ||  param.dutyRatio < 1
    eta = [];

    for i = 1: N

        [R,V] = GEq2RV(xg(i,1:6),1e-13, param.mu,param.J2, param.Re);

        % convert xguess from cartesian to coe.
        coe = CoordConv.vec2orbElem(R,V,param.mu);

       % fix this to represent the real time.
        theta_centre = wrapTo2Pi(eclipse_centre(param.tvec(i),coe,param));

        theta_centre_opp = wrapTo2Pi(theta_centre+pi);

        ArgofLat = wrapTo2Pi(coe(5) + coe(6));

        q1 = acos(cos(ArgofLat-theta_centre));
        q2 = acos(cos(ArgofLat-theta_centre_opp));

        if q1 < pi/4 || q2 < pi/4
            eta(i) = 0;
        else
            eta(i) = 1;
        end
    end
else
    eta = ones(1,N);


end
 % eta = ones(1,N);

if m1 == 1 %set as low thrust max lim
    % find dv limit by integrating Tmax/M using trapezium rule.
    dvmaglim =[];
    for i = 1:N
        a = param.Tmax/mass(i);
        b = param.Tmax/mass(i+1);
        h = param.dt;
        area = (a+b)/2*h;
        dvmaglim(i) =  param.T/param.m0*param.dt+ param.adddv;

    end


    dvmagarray =[];
    for i = 1: N
        dvmagarray = [dvmagarray, dvmaglim(i),dvmaglim(i),dvmaglim(i)];
    end
    %
    prob.blx = [ prob.blx, -dvmagarray,  zeros(1, N)];
    prob.bux = [prob.bux ,  dvmagarray,  eta.*dvmaglim];


elseif m1 == 2 %% get dv from initial guess

    dvmag = [];
    dvcomp = [];
    for i = 1: N

        dvmagnitude =  (norm(xg(i,7:9))+param.adddv);
        dvmag = [dvmag, eta(i)*dvmagnitude];
        dvcomp = [dvcomp, dvmagnitude*[1,1,1]];


    end

    prob.bux  = [prob.bux ,dvcomp, dvmag];
    prob.blx = [prob.blx,-dvcomp,zeros(1,N)];


else % set limits to inf

    dvmag = [];
    dvcomp = [];
    for i = 1: N

        dvmagnitude = 10; 
        dvmag = [dvmag, dvmagnitude];
        dvcomp = [dvcomp,  eta(i)*dvmagnitude*[1,1,1]];


    end

    prob.bux  = [prob.bux ,dvcomp, dvmag];
    prob.blx = [prob.blx,-dvcomp,zeros(1,N)];


%     prob.blx = [ prob.blx, -inf*ones(1,3*N),  zeros(1,N)];
%     prob.bux = [prob.bux ,  inf*ones(1,3*N),  inf*ones(1,N)];

end

%% 

prob.cones.type = zeros(N,1);

l = 1:3:3*N;
m = 2:3:3*N;
n = 3:3:3*N;

prob.cones.sub = [];
index = 6*N+6;
for i = 1:N % For each node
  
    prob.cones.sub  = [prob.cones.sub 9*N+6+i index+l(i) index+m(i) index+n(i)]; % Add the cone constraint
end
prob.cones.subptr = 1:4:4*(N); % Specify beginning index of each cone

[r,res]=mosekopt('minimize',prob,paramMosek);

% Display the primal solution.

x= res.sol.itr.xx';

if strcmp(res.sol.itr.solsta, 'PRIMAL_INFEASIBLE_CER') ||  strcmp(res.sol.itr.solsta, 'PRIMAL_ILLPOSED_CER')  %||  strcmp(res.sol.itr.solsta, 'UNKNOWN')
   changenu = true;
    'here'; 
else
    changenu = false;
end

end