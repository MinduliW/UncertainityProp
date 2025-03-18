function f = gen_kepler(K,p1,p2,varL)
% Generalized Kepler's equation. Appears as Eq. (26) in Ref:
%   Baù, G., Hernando-Ayuso, J. and Bombardelli, C., 2021.
%   "A generalization of the equinoctial orbital elements" 
%   Celestial Mechanics and Dynamical Astronomy, 133(11), pp.1-29.
%--------------------------------------------------------------------------
%
f = K + p1*cos(K) - p2*sin(K) - varL;
end

