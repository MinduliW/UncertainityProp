
function [A,B,C, D,dBdx] = getAandBMatrices(x, params)

% gravitational parameter
mu = params.mu;

p = x(1); f = x(2);g = x(3); 
h = x(4); k = x(5); L = x(6);

sinL = sin(L);
cosL = cos(L);
sqrtp_mu = sqrt(p/mu);

sqrt_mu = sqrt(mu);
q=1+f*cosL+g*sinL; 


s2 = 1 + h^2 + k^2; 

% A matrix
A = [ 0 0 0 0 0 sqrt(mu*p)*(q/p)^2]';

% B matrix
B = [0                      2*p/q*sqrtp_mu                0;
    sqrtp_mu*sinL  sqrtp_mu/q*((q+1)*cosL+f) -sqrtp_mu*g/q*(h*sinL-k*cosL);
    -sqrtp_mu*cosL sqrtp_mu*1/q*((q+1)*sinL+g) sqrtp_mu*f/q*(h*sinL-k*cosL);
    0                      0                                   sqrtp_mu*s2*cosL/(2*q);
    0                      0                                   sqrtp_mu*s2*sinL/(2*q);
    0                      0                                   sqrtp_mu*(h*sinL-k*cosL)/q];

% C matrix
C = [zeros(5,6); sqrt_mu*q^2*(-3/2)*p^(-5/2),...
    2*sqrt_mu*(1+ f*cosL + g*sinL)*cosL*p^(-3/2),...
    2*sqrt_mu*(1+ f*cosL + g*sinL)*sinL*p^(-3/2),...
    0,0, 2*sqrt_mu*(1+ f*cosL + g*sinL)*(-f*sinL + g*cosL)*p^(-3/2)];
 
% D matrix
D = A - C*x(1:6);

dBdx = getDBdx(x, params);


function dBdx = getDBdx(x, params)

mu = params.mu;

p = x(1); f = x(2);g = x(3); 
h = x(4); k = x(5); L = x(6);

sinL = sin(L);
cosL = cos(L);
sqrtp_mu = sqrt(p/mu);


dBdx(:,:,1) =                      [    0, (2*sqrtp_mu)/(f*cosL + g*sinL + 1) + p/(mu*sqrtp_mu*(f*cosL + g*sinL + 1)),                                                                          0;
 sinL/(2*mu*sqrtp_mu),       (f + cosL*(f*cosL + g*sinL + 2))/(2*mu*sqrtp_mu*(f*cosL + g*sinL + 1)),    (g*(k*cosL - h*sinL))/(2*mu*sqrtp_mu*(f*cosL + g*sinL + 1));
-cosL/(2*mu*sqrtp_mu),       (g + sinL*(f*cosL + g*sinL + 2))/(2*mu*sqrtp_mu*(f*cosL + g*sinL + 1)),   -(f*(k*cosL - h*sinL))/(2*mu*sqrtp_mu*(f*cosL + g*sinL + 1));
                         0,                                                                                          0, (cosL*(h^2 + k^2 + 1))/(2*mu*sqrtp_mu*(2*f*cosL + 2*g*sinL + 2));
                         0,                                                                                          0, (sinL*(h^2 + k^2 + 1))/(2*mu*sqrtp_mu*(2*f*cosL + 2*g*sinL + 2));
                        0,                                                                                          0,       -(k*cosL - h*sinL)/(2*mu*sqrtp_mu*(f*cosL + g*sinL + 1))];
 
 
dBdx(:,:,2) = [0,                                                                                             -(2*p*cosL*sqrtp_mu)/(f*cosL + g*sinL + 1)^2,                                                                                                                                          0;
0, ((cosL^2 + 1)*sqrtp_mu)/(f*cosL + g*sinL + 1) - (cosL*(f + cosL*(f*cosL + g*sinL + 2))*sqrtp_mu)/(f*cosL + g*sinL + 1)^2,                                                                 -(g*cosL*(k*cosL - h*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1)^2;
0,  (cosL*sinL*sqrtp_mu)/(f*cosL + g*sinL + 1) - (cosL*(g + sinL*(f*cosL + g*sinL + 2))*sqrtp_mu)/(f*cosL + g*sinL + 1)^2, (f*cosL*(k*cosL - h*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1)^2 - ((k*cosL - h*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1);
0,                                                                                                                                                  0,                                                                 -(2*cosL^2*sqrtp_mu*(h^2 + k^2 + 1))/(2*f*cosL + 2*g*sinL + 2)^2;
0,                                                                                                                                                  0,                                                            -(2*cosL*sinL*sqrtp_mu*(h^2 + k^2 + 1))/(2*f*cosL + 2*g*sinL + 2)^2;
0,                                                                                                                                                  0,                                                                    (cosL*(k*cosL - h*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1)^2];
 
 
dBdx(:,:,3) = [0,                                                                                             -(2*p*sinL*sqrtp_mu)/(f*cosL + g*sinL + 1)^2,                                                                                                                                          0;
0,  (cosL*sinL*sqrtp_mu)/(f*cosL + g*sinL + 1) - (sinL*(f + cosL*(f*cosL + g*sinL + 2))*sqrtp_mu)/(f*cosL + g*sinL + 1)^2, ((k*cosL - h*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1) - (g*sinL*(k*cosL - h*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1)^2;
0, ((sinL^2 + 1)*sqrtp_mu)/(f*cosL + g*sinL + 1) - (sinL*(g + sinL*(f*cosL + g*sinL + 2))*sqrtp_mu)/(f*cosL + g*sinL + 1)^2,                                                                  (f*sinL*(k*cosL - h*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1)^2;
0,                                                                                                                                                  0,                                                            -(2*cosL*sinL*sqrtp_mu*(h^2 + k^2 + 1))/(2*f*cosL + 2*g*sinL + 2)^2;
0,                                                                                                                                                  0,                                                                 -(2*sinL^2*sqrtp_mu*(h^2 + k^2 + 1))/(2*f*cosL + 2*g*sinL + 2)^2;
0,                                                                                                                                                  0,                                                                    (sinL*(k*cosL - h*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1)^2];
 
 
dBdx(:,:,4) =[0, 0,0;
0, 0,      -(g*sinL*sqrtp_mu)/(f*cosL + g*sinL + 1);
0, 0,       (f*sinL*sqrtp_mu)/(f*cosL + g*sinL + 1);
0, 0, (2*h*cosL*sqrtp_mu)/(2*f*cosL + 2*g*sinL + 2);
0, 0, (2*h*sinL*sqrtp_mu)/(2*f*cosL + 2*g*sinL + 2);
0, 0,  (sinL*sqrtp_mu)/(f*cosL + g*sinL + 1)];
 
 
dBdx(:,:,5) = [0, 0, 0;
0, 0,(g*cosL*sqrtp_mu)/(f*cosL + g*sinL + 1);
0, 0, -(f*cosL*sqrtp_mu)/(f*cosL + g*sinL + 1);
0, 0, (2*k*cosL*sqrtp_mu)/(2*f*cosL + 2*g*sinL + 2);
0, 0, (2*k*sinL*sqrtp_mu)/(2*f*cosL + 2*g*sinL + 2);
0, 0,        -(cosL*sqrtp_mu)/(f*cosL + g*sinL + 1)];
 
 
dBdx(:,:,6) =[0,                                                                                                                                                -(2*p*(g*cosL - f*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1)^2,                                                                                                                                                                       0;
cosL*sqrtp_mu, ((cosL*(g*cosL - f*sinL) - sinL*(f*cosL + g*sinL + 2))*sqrtp_mu)/(f*cosL + g*sinL + 1) - ((f + cosL*(f*cosL + g*sinL + 2))*(g*cosL - f*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1)^2,           - (g*(h*cosL + k*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1) - (g*(g*cosL - f*sinL)*(k*cosL - h*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1)^2;
sinL*sqrtp_mu, ((cosL*(f*cosL + g*sinL + 2) + sinL*(g*cosL - f*sinL))*sqrtp_mu)/(f*cosL + g*sinL + 1) - ((g + sinL*(f*cosL + g*sinL + 2))*(g*cosL - f*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1)^2,             (f*(h*cosL + k*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1) + (f*(g*cosL - f*sinL)*(k*cosL - h*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1)^2;
                  0,                                                                                                                                                                                                                    0, - (sinL*sqrtp_mu*(h^2 + k^2 + 1))/(2*f*cosL + 2*g*sinL + 2) - (cosL*(2*g*cosL - 2*f*sinL)*sqrtp_mu*(h^2 + k^2 + 1))/(2*f*cosL + 2*g*sinL + 2)^2;
                  0,                                                                                                                                                                                                                    0,   (cosL*sqrtp_mu*(h^2 + k^2 + 1))/(2*f*cosL + 2*g*sinL + 2) - (sinL*(2*g*cosL - 2*f*sinL)*sqrtp_mu*(h^2 + k^2 + 1))/(2*f*cosL + 2*g*sinL + 2)^2;
                  0,                                                                                                                                                                                                                    0,                 ((h*cosL + k*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1) + ((g*cosL - f*sinL)*(k*cosL - h*sinL)*sqrtp_mu)/(f*cosL + g*sinL + 1)^2];
 

