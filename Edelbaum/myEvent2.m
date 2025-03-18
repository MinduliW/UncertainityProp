
function [value, isterminal, direction] = myEvent2(time,states,RAANexpt,param)

% convert states to coe.
coe = CoordConv.mee2coe(states);
hill = kep2hill(coe, param.mu);
hillMean = osculating2meanHill(hill, param.mu, param.J2, param.Re);
mean_waiting = hill2kep(hillMean, param.mu);


value      = abs(mean_waiting(4) - RAANexpt)<0.1*pi/180;
isterminal = 1;   % Stop the integration
direction  = 0;
end

