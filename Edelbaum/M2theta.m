function  theta = M2theta(M,e)
E = M;
if e<1
for i=1:20;
    ddf = (e*cos(E)-1);
    E  = E-(M-E + e*sin(E))/ddf;
    theta = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
end

else 
for i=1:20;
    ddf = (1 - e*cosh(E));
    E  = E-(M + E - e*sinh(E))/ddf;
    theta  = 2*atan(sqrt((1+e)/(e-1))*tanh(E/2));
end
end