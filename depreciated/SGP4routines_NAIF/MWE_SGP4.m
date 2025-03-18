clear;
clc;

constants;

addpath('/Users/minduli/libraries_mice/mice/src/mice/')
addpath('/Users/minduli/libraries_mice/mice/lib')
addpath('../parametric analysis')
cspice_furnsh('/Users/minduli/Astroscale_ADR/Main/SGP4routines_NAIF/kernel.txt')

Norad1 = '33500';
Norad2 = '39766';
Norad3 = '33492';

SatID1 = num2str(-100000-str2double(Norad1));
SatID2 = num2str(-100000-str2double(Norad2));
SatID3 = num2str(-100000-str2double(Norad3));


datetimestring = '2022-03-25 06:37:13';
t0 = cspice_str2et(datetimestring);
jd1s = cspice_et2utc(t0,'J', 6);

newStr = split(jd1s,{' '});
jd1 = str2double(newStr{2});

%%
%To print to screen the window of validity of kernels (example on Norad1):
SPK1 = strcat('SPK',Norad1,'.bsp');
MAXIV  = 1000;
WINSIZ = 2 * MAXIV;
ids = cspice_spkobj( SPK1, MAXIV );
cover = cspice_spkcov( SPK1, ids(1), WINSIZ );

sprintf('Initial epoch: %s', cspice_et2utc(cover(1),'C',2))
sprintf('Final epoch: %s', cspice_et2utc(cover(2),'C',2))

% timevec
timevec = jd1:0.1:jd1+1; 

RAANs = zeros(size(timevec));
RAAN_cur = zeros(size(timevec));
for i = 1: length(timevec)
    
    newstr = strcat('JD', {' '}, num2str(timevec(i)));
    
    datestring = cspice_str2et(newstr);
 
    State1 = cspice_spkezr(SatID1, datestring, 'J2000', 'none', '399'); 
    
    % convert to keplerian elements.
    x = CoordConv.vec2orbElem(State1(1:3),State1(4:6),3.986005e5);
    
    [rr, vv, E,epoch, mass] = getPosition(timevec(i), 'H2AF15', param.mu,...
        param.J2, param.Re);

    RAAN_cur(i) = E(4);
    RAANs(i) = x(5); 
    
end

figure; hold on;
plot(timevec-timevec(1), unwrap(RAANs)*180/pi);
p = plot(timevec-timevec(1), unwrap(RAAN_cur)*180/pi);
plot_latex(p, 'time (JD)', 'RAAN (deg)','' ,'', {'from SGP4', 'current implementation'})


%To retrieve state at specific epoch:
State1 = cspice_spkezr(SatID1, t0, 'J2000', 'none', '399'); %retrieve Sat1 at t0 in coordinates j2000 with no light correction centered at the earth (399)
x = CoordConv.vec2orbElem(State1(1:3),State1(4:6),3.986005e5);
    






