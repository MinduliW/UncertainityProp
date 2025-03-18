function density = Density_nrlmsise00(altitude, timeinJD, param)

longitude =rand(1)*180;
latitude =rand(1)*180-90;

datevec = datetime(timeinJD, 'ConvertFrom','juliandate', 'Format','dd-MM-yyy');

y = year(datevec);
dayOfYear = day(datevec,'dayofyear');
UTseconds =  second(datevec);

[ f107Average, f107Daily, magneticIndex ] = computeSW( param.SWmatDaily, param.SWmatMonthlyPred, timeinJD );

[~ ,density] = atmosnrlmsise00(altitude,latitude,longitude,y,dayOfYear,UTseconds,...
    f107Average, f107Daily, magneticIndex, 'Oxygen', 'None');
density = density(6);

end
