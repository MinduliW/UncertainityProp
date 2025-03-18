function hill = kep2hill(kep, mu)
if kep(2) == 0 
    kep(2) = 1e-5;
end


p = kep(1)*(1.0 - kep(2)*kep(2));
f = kep(6);
hill(5) = sqrt(mu*p);
hill(1) = p/(1.0 + kep(2)*cos(f));
hill(2) = f + kep(5);
hill(3) = kep(4);
hill(4) = (hill(5)/p)*kep(2)*sin(f);
hill(6) = hill(5)*cos(kep(3));

end





