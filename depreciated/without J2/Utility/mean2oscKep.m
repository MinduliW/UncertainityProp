function x0KepOsc = mean2oscKep(x0Kep, mu, J2, Re)

hill = kep2hill(x0Kep, mu);
hillosc = mean2osculatingHill(hill, mu, J2, Re);
x0KepOsc = hill2kep(hillosc,mu);

end
