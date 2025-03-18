function EA = true2eccAnomaly(theta , e)
EA = 2.0 * atan2(sqrt(1. - e)*sin(theta / 2.), sqrt(1. + e) * cos(theta / 2.));

end

