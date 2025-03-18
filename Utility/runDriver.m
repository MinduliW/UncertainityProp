function [As,Bs, constants, TR, xfconv,nlIndx] = runDriver(paramArray, xguess)

[STMs, TR, xfconversion,nlIndx] = driver(paramArray, xguess);
As = STMs(:,1:6);
Bs = STMs(:,7:9);
constants = STMs(:,10);
xfconv.A = xfconversion(:,1:6);
xfconv.B = xfconversion(:,7);
end
