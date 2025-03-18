function [x0, ub,lb] = getInitGuess(debrisID, param)
 
% Any date works here, it's just to get a and inc.
datevec = '19-June-2024';
t0 = juliandate(datevec,'dd-mmm-yyyy'); %+ 1000;


debris = zeros(length(debrisID),6);

for i = 1:length(debrisID)
    [~, ~, debris(i,:)] = getPosition(t0 ,debrisID(i),param.mu, param.J2, param.Re);
end

% Do not allow waiting orbits to have an altitude below the one of the rocketlab
a_RL = param.target_altitude + param.Re;
v_RL = sqrt(param.mu/a_RL);

x0 = [];
lb = [];
ub = []; 

for j = 2: length(debrisID)
    a0 = param.target_altitude + param.Re; inc0 = debris(j-1,3);
    v0 = sqrt(param.mu/a0);
    
   af =  debris(j,1);incf = debris(j,3);
   vf = sqrt(param.mu/af);
   
   
x0 = [x0, v0, inc0];
lb = [lb , min([v0,vf])-250, min([inc0,incf])-1*180/pi];
ub = [ub, v_RL+250, max([inc0,incf])+1*180/pi];

end


%% add time limits
x0 = [x0, t0];
lb = [lb, t0-500];
ub = [ub, t0+500];
