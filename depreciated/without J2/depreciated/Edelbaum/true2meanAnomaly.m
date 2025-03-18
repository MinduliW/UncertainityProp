function MA = true2meanAnomaly(theta, e)

	E = true2eccAnomaly(theta, e);

	MA = E - e*sin(E);
end
