
function [fsolveoptions, odeoptions,bvpoptions] = fsolOpt(tol)
fsolveoptions = optimoptions('fsolve','Display','iter');
fsolveoptions.FunctionTolerance = tol;
fsolveoptions.FunctionTolerance = tol;
fsolveoptions.MaxFunctionEvaluations = 6e3;
odeoptions =  odeset('RelTol', tol,'AbsTol',tol);
bvpoptions = bvpset('RelTol', tol,'AbsTol',tol);

end
