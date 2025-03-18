function coemean = osc2meankep(coe, mu, J2, Re)

hill = kep2hill(coe, mu);
hillmean = osculating2meanHill(hill,mu, J2, Re);
coemean= hill2kep(hillmean, mu);


end
