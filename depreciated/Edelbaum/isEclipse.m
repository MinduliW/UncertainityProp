function eta = isEclipse(time, semiMajor, eccentricity, inclination, RAAN, AOP, theta, param)

t = time*param.TU/86400 + param.t0;
kepElements = [semiMajor, eccentricity, inclination, RAAN, AOP,theta];
rr = CoordConv.po2pv(kepElements, param.mu); %m, m/s
r_Earth2Sun = getSunPosVec(t, param);

eta = eclipse(-r_Earth2Sun+rr, rr(1:3), param);

