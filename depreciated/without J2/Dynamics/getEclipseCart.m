function v = getEclipseCart(t, x, param)
%Shadow function computation
[~,r_Earth2Sun] = getSunPosVec(t + 2400000.5+60000+param.deltat0,param);

%r_Earth2Sun = sunPosWang(t,param);

posandvel = x(1:6); % CoordConv.ep2pv(x(1:6),param.mu);
rEarthToSat = posandvel(1:3);

rSat2Earth = -rEarthToSat;
rSat2Sun = rSat2Earth+r_Earth2Sun;

%v =  eclipse(-rSat2Sun, rEarthToSat,param);
v= 1- Smoothed_eclipse(param.cs,param.ct, rSat2Sun ,...
    rSat2Earth,param);
end
