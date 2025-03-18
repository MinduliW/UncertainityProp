 function [coemean,coe,timevec] = errorcalc(times, states, af,incf,targetRAAN,param)


% coordinate conversion
coemean =[]; coe = []; 

% while isnan(states(end,1))
%     states = states(1:end-1,:);
%     times = times(1:end-1);
% end
% 
% 
% % if length(times)>  1e4
% %     times = times(1:100:end);
% %     states = states(1:100:end,:);
% % else
% %     step = 1; 
% % end

timevec = [];

for j = 1: 1: length(times)
    
    timevec =  [timevec ; times(j)];
    
    
    coec = CoordConv.mee2coe(states(j,1:6));
    coe = [coe; coec];
    
    % convert to hill
    hill = kep2hill(coec, param.mu);
    
    % calculate mean hill
    hillMean = osculating2meanHill(hill, param.mu, param.J2, param.Re);
    
    cmean = hill2kep(hillMean, param.mu);
    
    coemean = [coemean; cmean];
end




% Accuracy
meansmadiff= (coemean(end,1)-af);
meanincdiff = (coemean(end,3)-incf)*180/pi; 
meanOmegadiff = wrapTo2Pi(unwrap(wrapTo2Pi(coemean(end,4)))-wrapTo2Pi(targetRAAN))*180/pi; 
meanediff = coemean(end,2);


fprintf('Error in mean sma: %f km \n', meansmadiff*param.LU/1e3);
fprintf('Error in mean inc: %f deg \n', meanincdiff);
fprintf('Error in mean RAAN: %f deg \n', meanOmegadiff);
fprintf('Error in mean ecc: %f  \n', meanediff);

end
