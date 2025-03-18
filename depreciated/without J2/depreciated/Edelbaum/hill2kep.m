function kep = hill2kep(hill, mu)

r  = hill(1);
th = hill(2);
nu = hill(3);
R  = hill(4);
Th = hill(5);
Nu = hill(6);
    
i = acos(Nu/Th);

if ~isreal(i)
    'here'
    i = pi/2;
    nu = 0;
end

cs  =   (-1.0 + pow(Th,2)/(mu*r))*cos(th) + (R*Th*sin(th))/mu;
ss  =  -((R*Th*cos(th))/mu) + (-1.0 + pow(Th,2)/(mu*r))*sin(th);
e = sqrt(cs*cs+ss*ss);


p = Th*Th/mu;
if e == 0 
    e = 1e-5;
end


costrue = 1/e*(p/r-1.0);

if abs(costrue) >1
    costrue = sign(costrue)*1; 
end


f = acos(costrue);
    

if ((R)<0.0)
    f = 2.0*pi-f;
end

a = p/(1-e*e);

kep(1) = a;
kep(2) = e;
kep(3) = wrapTo2Pi(i);
kep(4) = wrapTo2Pi(nu);
kep(5) = wrapTo2Pi(th-f);
kep(6) = wrapTo2Pi(f);


end
