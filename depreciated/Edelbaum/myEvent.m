
function [value, isterminal, direction] = myEvent(time,states,smaexpt,param)

coe = CoordConv.mee2coe(states);
hill = kep2hill(coe, param.mu);
hillMean = osculating2meanHill(hill, param.mu, param.J2, param.Re);
mean_waiting = hill2kep(hillMean, param.mu);

value      = abs(mean_waiting(1) - smaexpt)<0.5e3;
isterminal = 1;   % Stop the integration
direction  = 0;
end


