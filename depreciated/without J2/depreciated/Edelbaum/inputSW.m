function [ SWmatDaily, SWmatMonthlyPred ] = inputSW( swfName )
%READSW reads space weather file from CelesTrack
% [  ] = READSW(SWFNAME, JDATE, UTHR)
%
% Inputs for INPUTSW are:
% SWFNAME   :a string that contains space weather name
%
% Outputs for INPUTSW are:
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
%

%% Internal variables definition
% line length
linelenght = 132;

% Julian Date of 1957 10 01 00:00:00
jdate1957 = 2436112.5;

% Julian Date of 2000 01 01 00:00:00
jdate2000 =  2451544.5;

%% File processing

fid = fopen(swfName,'r+');

% Get rid of initial lines
for i=1:15
    fgetl(fid);
end

% Read number of observed points
str = fgetl(fid);
obs_pnt = str2double(str(21:25));

% get rid of BEGIN OBSERVED
fgetl(fid);

% Jump to line corresponding to 2000 01 01 00:00:00
dataline = floor(jdate2000-jdate1957)*linelenght;
fseek(fid, dataline, 0);

n_daily_obs = obs_pnt-floor(jdate2000-jdate1957);

SWaux = zeros(n_daily_obs, 11);

for i = 1:n_daily_obs
    
    str = fgetl(fid);
    
    SWaux(i, 1) = str2double(str(94:98)); % F10.7 Daily
    
    SWaux(i, 2) = str2double(str(102:106)); % F10.7 Average
    
    SWaux(i, 3) = str2double(str(80:82)); % Daily Magnetic index
    
    SWaux(i, 4:11) = str2num([str(47:50),str(51:54),str(55:58),str(59:62),...
                         str(63:66),str(67:70),str(71:74),str(75:78)]); % Daily 3h APs
                     
end

for i=1:3
    str=fgetl(fid);
end

pdt_pnt = str2double(str(28:29));

SWmatDaily = zeros( n_daily_obs + pdt_pnt, 11);
SWmatDaily(1:n_daily_obs, :) = SWaux;

clear SWaux;

% get rid of BEGIN DAILY_PREDICTED
fgetl(fid);

for i = n_daily_obs+1:n_daily_obs+pdt_pnt
    
    str = fgetl(fid);
    
    SWmatDaily(i, 1) =  str2double(str(94:98)); % F10.7 Daily
    
    SWmatDaily(i, 2) = str2double(str(102:106)); % F10.7 Average
    
    SWmatDaily(i, 3) = str2double(str(80:82)); % Daily Magnetic index
    
    SWmatDaily(i, 4:11) = str2num([str(47:50),str(51:54),str(55:58),str(59:62),...
                         str(63:66),str(67:70),str(71:74),str(75:78)]); % Daily 3h APs
    
end

for i=1:3
    str=fgetl(fid);
end

mpd_pnt = str2double(str(30:31));

SWmatMonthlyPred = zeros(mpd_pnt, 2);

% get rid of BEGIN MONTHLY_PREDICTED
fgetl(fid);

for i=1:mpd_pnt
    str=fgetl(fid);
    
    SWmatMonthlyPred(i, 1) =  str2double(str(94:98)); % F10.7 Daily
    
    SWmatMonthlyPred(i, 2) = str2double(str(102:106)); % F10.7 Average
    
    %SWmatMonthlyPred(i, 3:11) = 4*ones(1,9); % Daily Magnetic index

end

fclose(fid);

end