function meanCart = osc2meanCart(rr,vv,mu, J2, Re)

osckep = CoordConv.vec2orbElem(rr,vv,mu);

meankep = osc2meankep(osckep, mu, J2, Re);

[rr,vv] = CoordConv.po2pv(meankep, mu);
meanCart = [rr;vv];
