function  theta = theta2M(theta,e)


if e > 1

H = 2.0 * atan2(sqrt(e -1)*sin(theta / 2.), sqrt(1. + e) * cos(theta / 2.));
theta = H - e*sinh(H);

else

E = 2.0 * atan2(sqrt(1. - e)*sin(theta / 2.), sqrt(1. + e) * cos(theta / 2.));
theta = E - e*sin(E);
end

end