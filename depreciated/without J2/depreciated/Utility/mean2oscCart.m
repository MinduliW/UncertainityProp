function oscCart = mean2oscCart(rr,vv,mu, J2, Re)

meankep = CoordConv.vec2orbElem(rr,vv,mu);
osckep = mean2oscKep(meankep, mu, J2, Re);

[rr,vv] = CoordConv.po2pv(osckep, mu);
oscCart = [rr;vv];
