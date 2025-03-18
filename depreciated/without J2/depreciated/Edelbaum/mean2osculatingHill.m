

function hillOsc = mean2osculatingHill(hillMean, param)

    %Mean Delaunay elements
    r  = hillMean(1); 
    th = hillMean(2); 
    nu = hillMean(3);
    R  = hillMean(4); 
    Th = hillMean(5); 
    Nu = hillMean(6); 
    
    ci = Nu/Th;
    si = sqrt(1.0-ci*ci);
    cs  =   (-1.0 + pow(Th,2)/(param.mu*r))*cos(th) + (R*Th*sin(th))/param.mu;
    ss  =  -((R*Th*cos(th))/param.mu) + (-1.0 + pow(Th,2)/(param.mu*r))*sin(th);
    e = sqrt(cs*cs+ss*ss);
    eta  = sqrt(1.0-e*e);
    beta = 1.0/(1.0+eta);
    
    p = Th*Th/param.mu;
    costrue = 1/e*(p/r-1);
    
    f = acos(costrue);
    
    if (R<0.0) 
        f = 2.0*pi-f;
    end
    
    
    M = true2meanAnomaly(f,e);
    
    phi  = f - M;
    
    r = r -  (1.0)*((pow(param.Re,2)*beta*param.J2)/(2.*r) - (3*pow(param.Re,2)*beta*param.J2*pow(si,2))/(4.*r) +...
                     (pow(param.Re,2)*eta*param.J2*pow(param.mu,2)*r)/pow(Th,4) - (3*pow(param.Re,2)*eta*param.J2*pow(param.mu,2)*r*pow(si,2))/(2.*pow(Th,4)) +...
                     (pow(param.Re,2)*param.J2*param.mu)/(2.*pow(Th,2)) - (pow(param.Re,2)*beta*param.J2*param.mu)/(2.*pow(Th,2)) -...
                     (3.*pow(param.Re,2)*param.J2*param.mu*pow(si,2))/(4.*pow(Th,2)) + (3*pow(param.Re,2)*beta*param.J2*param.mu*pow(si,2))/(4.*pow(Th,2)) -...
                     (pow(param.Re,2)*param.J2*param.mu*pow(si,2)*cos(2*th))/(4.*pow(Th,2)));
    
    
    th = th - (1.0)*((-3.*pow(param.Re,2)*param.J2*pow(param.mu,2)*phi)/pow(Th,4) + (15.*pow(param.Re,2)*param.J2*pow(param.mu,2)*phi*pow(si,2))/(4.*pow(Th,4)) -...
                      (5.*pow(param.Re,2)*param.J2*param.mu*R)/(2.*pow(Th,3)) - (pow(param.Re,2)*beta*param.J2*param.mu*R)/(2.*pow(Th,3)) +...
                      (3.*pow(param.Re,2)*param.J2*param.mu*R*pow(si,2))/pow(Th,3) + (3.*pow(param.Re,2)*beta*param.J2*param.mu*R*pow(si,2))/(4.*pow(Th,3)) -...
                      (pow(param.Re,2)*beta*param.J2*R)/(2.*r*Th) + (3.*pow(param.Re,2)*beta*param.J2*R*pow(si,2))/(4.*r*Th) +...
                      (-(pow(param.Re,2)*param.J2*param.mu*R)/(2.*pow(Th,3)) + (pow(param.Re,2)*param.J2*param.mu*R*pow(si,2))/pow(Th,3))*cos(2.*th) +...
                      (-(pow(param.Re,2)*param.J2*pow(param.mu,2))/(4.*pow(Th,4)) + (5.*pow(param.Re,2)*param.J2*pow(param.mu,2)*pow(si,2))/(8.*pow(Th,4)) +...
                       (pow(param.Re,2)*param.J2*param.mu)/(r*pow(Th,2)) - (3.*pow(param.Re,2)*param.J2*param.mu*pow(si,2))/(2.*r*pow(Th,2)))*sin(2.*th));
    
    nu = nu - (1.0)*((3.*pow(param.Re,2)*ci*param.J2*pow(param.mu,2)*phi)/(2.*pow(Th,4)) + (3.*pow(param.Re,2)*ci*param.J2*param.mu*R)/(2.*pow(Th,3)) +...
                      (pow(param.Re,2)*ci*param.J2*param.mu*R*cos(2.*th))/(2.*pow(Th,3)) +...
                      ((pow(param.Re,2)*ci*param.J2*pow(param.mu,2))/(4.*pow(Th,4)) - (pow(param.Re,2)*ci*param.J2*param.mu)/(r*pow(Th,2)))*sin(2.*th));
    
    
    R = R  - (1.0)*(-(pow(param.Re,2)*beta*param.J2*R)/(2.*pow(r,2)) + (3.*pow(param.Re,2)*beta*param.J2*R*pow(si,2))/(4.*pow(r,2)) -...
                     (pow(param.Re,2)*eta*param.J2*pow(param.mu,2)*R)/(2.*pow(Th,4)) + (3.*pow(param.Re,2)*eta*param.J2*pow(param.mu,2)*R*pow(si,2))/(4.*pow(Th,4)) +...
                     (pow(param.Re,2)*param.J2*param.mu*pow(si,2)*sin(2.*th))/(2.*pow(r,2)*Th));
    
    
    Th = Th  - (1.0)*(((pow(param.Re,2)*param.J2*pow(param.mu,2)*pow(si,2))/(4.*pow(Th,3)) - (pow(param.Re,2)*param.J2*param.mu*pow(si,2))/(r*Th))*cos(2.*th) -...
                       (pow(param.Re,2)*param.J2*param.mu*R*pow(si,2)*sin(2.*th))/(2.*pow(Th,2)));
    
    Nu = Nu +  0;
    
    
    hillOsc(1) = r;
    hillOsc(2) = th;
    hillOsc(3) = nu;
    hillOsc(4) = R;
    hillOsc(5) = Th;
    hillOsc(6) = Nu;
end

