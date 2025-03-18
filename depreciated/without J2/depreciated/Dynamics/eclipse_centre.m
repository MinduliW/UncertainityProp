function centre = eclipse_centre(t,coe,param)

%centre = 0;
norbits = 2; 
period =  2*pi*sqrt(coe(1)^3/param.mu);

TA = wrapTo2Pi(coe(end));

% on a grid of theta, determine the start and end (hence centre) of the
% eclipse for the orbit.
theta = (linspace(TA,TA+norbits*2*pi,30));
time = linspace(t,t+norbits*period,30);

% for each theta, get the eclipse.
for i = 1: length(theta)
    ecl(i) = isEclipse(time(i), coe(1),coe(2), coe(3), coe(4), coe(5), theta(i), param);
end

% get eclipse thetas
eclipseindx = find(ecl ==0);


if isempty(eclipseindx) 
    centre = 0;
elseif length(eclipseindx)<2
    centre = theta(eclipseindx);
else
    % find the largest gap present.
    for i = 1: length(eclipseindx)-1
        
        diff(i) = eclipseindx(i+1) - eclipseindx(i);
        
    end
    
    relt = find(diff >1);
    
    if length(relt)>1
        centre = (theta(eclipseindx(relt(1)+1))+theta(eclipseindx(relt(2))))/2;
    elseif isempty(relt)
        
        % this happens if the orbit is fully in sunlight.
        centre = 0;
        
    else
        
        if length(eclipseindx)<1 || length(relt)<1 
            'error in eclipse: check eclipse_centre';
        end
        
            

        centre = (theta(eclipseindx(1))+theta(eclipseindx(relt(1))))/2;
    end
    
    centre  =(centre);
end

%argLatCen = centre +  wrapTo2Pi(coe(4));


% figure; hold on;
% plot((theta)*180/pi, ecl);
% plot(centre*ones(size(ecl))*180/pi, linspace(0,1, length(ecl)));

end