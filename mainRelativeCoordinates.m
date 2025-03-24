%% To do:
% Define satellite 1: Cheif , at [500 km , 0, 98 deg,0,0,0 ]
% Define satellite 2: Deputy , at [500 km , 0, 98 deg,0,0,10 deg]

% Convert both from Keplerian to Cartesian

% Propagate both in cartesian. 

% Calculate the relative LVLH state at the start and at the end. 


% Introduce errors to both the cheif and deputy like before.

% Convert perturbed states to LVLH and ROE. 

 %Propagate each. 

 % Check the uncertainity evolution. 

 % PS: I will give you the propagation and conversion code for everything
 % in Matlab, you might have to adopt them to Cpp to calculate the STMs to
 % calculate the nonlinearity index. 
 % Rel2Cart2.m and Cart2Rel2.m - from Cartesian to LVLH and back 
 % roe2kep and kep2roe - ROE to keplerian and back
 % dynamics_roe.m - propagation in ROE. 
 % NonlinearRel - nonlinear LVLH prop - also look at Clohessy Wlitshire
 % dynamics here. 

 % Sources 
 % Spacecraft formation flying textbook - covers all you need to know about
 % relative motion
 % Relative orbital elements -
 % https://arc.aiaa.org/doi/pdf/10.2514/1.G006175 - this is a helpful
 % source