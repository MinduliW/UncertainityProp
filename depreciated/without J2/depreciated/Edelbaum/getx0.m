function diff = getx0(theta, a0, RAAN0, inc0,param)

x0 = CoordConv.kepler2MEOE([a0, 0,0,RAAN0,inc0,theta]);

% get the mean elements.
coe = CoordConv.mee2coe(x0);


hill = kep2hill(coe, param.mu);
hillMean = osculating2meanHill(hill, param.mu, param.J2, param.Re);
coemean = hill2kep(hillMean, param.mu);



diff(1) = coemean(1)-a0; 
%diff(2) = coemean(3)- inc0;
%diff(3) = coemean(4) - RAAN0;


end
