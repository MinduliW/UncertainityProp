function hillMean = osculating2meanHill(hillOsc, mu, J2, rE)

% Mean Delaunay elements
r   = hillOsc(1);
th  = hillOsc(2);
nu  = hillOsc(3);
R   = hillOsc(4);
Th  = hillOsc(5);
Nu  = hillOsc(6);


ci = Nu/Th;
si = sqrt(1.0-ci*ci);

rEsquared = pow(rE,2);
musquared = pow(mu,2);
Thsquared = pow(Th,2);
sisquared = pow(si,2);
Thcubed = pow(Th,3);
Thfour = pow(Th,4);


cs  =   (-1.0 + Thsquared/(mu*r))*cos(th) + (R*Th*sin(th))/mu;
ss  =  -((R*Th*cos(th))/mu) + (-1.0 + Thsquared/(mu*r))*sin(th);
e = sqrt(cs*cs+ss*ss);

if e == 0 
    e = 1e-5;
end

eta  = sqrt(1.0-e*e);


beta = 1.0/(1.0+eta);
p = Th*Th/mu;
costrue = 1/e*(p/r-1);

if abs(costrue) >1
    costrue = sign(costrue)*1; 
end



f = acos(costrue);

if ((R)<0.0)
    f = 2.0*pi-f;
end


M = theta2M(f,e);

phi  = f - M;

r = r +  ((rEsquared*beta*J2)/(2.*r) - (3*rEsquared*beta*J2*sisquared)/(4.*r) +...
(rEsquared*eta*J2*musquared*r)/Thfour - (3*rEsquared*eta*J2*musquared*r*sisquared)/(2.*Thfour) +...
(rEsquared*J2*mu)/(2.*Thsquared) - (rEsquared*beta*J2*mu)/(2.*Thsquared) -...
(3.*rEsquared*J2*mu*sisquared)/(4.*Thsquared) + (3*rEsquared*beta*J2*mu*sisquared)/(4.*Thsquared) -...
(rEsquared*J2*mu*sisquared*cos(2*th))/(4.*Thsquared));

th = th + ((-3.*rEsquared*J2*musquared*phi)/Thfour + (15.*rEsquared*J2*musquared*phi*sisquared)/(4.*Thfour) -...
(5.*rEsquared*J2*mu*R)/(2.*Thcubed) - (rEsquared*beta*J2*mu*R)/(2.*Thcubed) +...
(3.*rEsquared*J2*mu*R*sisquared)/Thcubed + (3.*rEsquared*beta*J2*mu*R*sisquared)/(4.*Thcubed) -...
(rEsquared*beta*J2*R)/(2.*r*Th) + (3.*rEsquared*beta*J2*R*sisquared)/(4.*r*Th) +...
(-(rEsquared*J2*mu*R)/(2.*Thcubed) + (rEsquared*J2*mu*R*sisquared)/Thcubed)*cos(2.*th) +...
(-(rEsquared*J2*musquared)/(4.*Thfour) + (5.*rEsquared*J2*musquared*sisquared)/(8.*Thfour) +...
(rEsquared*J2*mu)/(r*Thsquared) - (3.*rEsquared*J2*mu*sisquared)/(2.*r*Thsquared))*sin(2.*th));


nu = nu + ((3.*rEsquared*ci*J2*musquared*phi)/(2.*Thfour) + (3.*rEsquared*ci*J2*mu*R)/(2.*Thcubed) +...
(rEsquared*ci*J2*mu*R*cos(2.*th))/(2.*Thcubed) +...
((rEsquared*ci*J2*musquared)/(4.*Thfour) - (rEsquared*ci*J2*mu)/(r*Thsquared))*sin(2.*th));


R = R  + (-(rEsquared*beta*J2*R)/(2.*pow(r,2)) + (3.*rEsquared*beta*J2*R*sisquared)/(4.*pow(r,2)) -...
(rEsquared*eta*J2*musquared*R)/(2.*Thfour) + (3.*rEsquared*eta*J2*musquared*R*sisquared)/(4.*Thfour) +...
(rEsquared*J2*mu*sisquared*sin(2.*th))/(2.*pow(r,2)*Th));


Th = Th  + (((rEsquared*J2*musquared*sisquared)/(4.*Thcubed) - (rEsquared*J2*mu*sisquared)/(r*Th))*cos(2.*th) -...
(rEsquared*J2*mu*R*sisquared*sin(2.*th))/(2.*Thsquared));

Nu = Nu +  0;

hillMean(1) = r;
hillMean(2) = th;
hillMean(3) = nu;
hillMean(4) = R;
hillMean(5) = Th;
hillMean(6) = Nu;


end



    

