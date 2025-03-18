function [ f107Average, f107Daily, magneticIndex ] = computeSW( SWmatDaily, SWmatMonthlyPred, jdate )
%READSW reads space weather file from CelesTrack
% [  ] = READSW(SWFNAME, JDATE, UTHR)
%
% Inputs for READSW are:
% SWMATDAILY : 
%              matrix for F10.7Daily, F10.7Average, magnetic index
%              Daily observed and predicted AP (8)
%              from 1 Jan 2000 to end of Daily predicted
%
% SWMATMONTHLYPRED : 
%              matrix for Monthly predicted F10.7Daily, F10.7Average
%              Magnetic index and AP (8) from 1 Jan 2000 to end of 
%              Daily predicted
%
% DATE      : Julian Date
%
%
%

%% Outputs initializations (sets default values for atmosnrlmsise00)
magneticIndex = 4*ones(1,7);
f107Daily = 150;
f107Average = 150;

%% Internal variables definition

% Julian Date of 2000 01 01 00:00:00
jdate2000 =  2451544.5;

%% Determine UT hour
[UTyr,UTmo,~,UThr,~,~] = invjday(jdate);

%% File processing
auxMI = zeros(1,32);

row = floor(jdate-jdate2000)+1;

if row <= size(SWmatDaily,1)
    
    f107Daily =  SWmatDaily(row-1, 1);
    
    f107Average = SWmatDaily(row, 2);
    
    magneticIndex(1) = SWmatDaily(row, 3);
    column = ceil((24-UThr)/3);
    
    for i=1:4
        
        auxMI(i*8-7:i*8) = SWmatDaily(row-(4-i), 4:11);

    end
        
    auxMI = fliplr(auxMI);
    
    magneticIndex(2:5) = auxMI(column:column+3);
    
    magneticIndex(6) = mean(auxMI(column+4:column+11));
    magneticIndex(7) = mean(auxMI(column+12:column+19));
    
else
    
    % Determine UTyr and UTmo of Monthly prediction beginning
    [UTyrMP, UTmoMP, ~, ~, ~, ~] = invjday(jdate2000 + size(SWmatDaily,1));
    
     dUTyr = UTyr-UTyrMP;
     dUTmo = UTmo-UTmoMP;
        
     dmon = dUTmo + dUTyr*12;
        
     if dmon<0
         dmon=0;
     end
       
     if dmon<size(SWmatMonthlyPred,1);
        
        f107Daily = SWmatMonthlyPred(dmon+1,1);
        f107Average = SWmatMonthlyPred(dmon+1,2);
         
     end
        
end


end

