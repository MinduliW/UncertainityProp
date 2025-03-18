function w = shadowArc(semiMajor, inclination, RAAN, tSeg,param)


period = 2*pi*sqrt(semiMajor^3/param.mu);

theta = linspace(0, 2*pi, 100);
time = linspace(tSeg, tSeg + period, 100);

eta = [];


% Get sun vector at each position through 1 orbit.
for i = 1:length(time)
    
  
    eta(i) = isEclipse(time(i), semiMajor,0, inclination, RAAN,0, theta(i), param);
 
    if   ~isreal(eta(i))
        disp(eta(i));
      
        disp('Warning, Shadow function returns compelex value!')
        
    end
    
end

eta_eclipse = eta(eta == 0);
eclipse_times = length(eta_eclipse)/length(eta);

w = 1  - eclipse_times;

end


