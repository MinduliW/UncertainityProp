datetimestring = '2022-03-25 06:37:13';

test1 = datetime(datetimestring);

juliandate(test1)


t0 = cspice_str2et(datetimestring);
jd1s = cspice_et2utc(t0,'J', 6);

newStr = split(jd1s,{' '});
jd1 = str2double(newStr{2});

%jd1 = jd1+1; 

newstr = strcat('JD', {' '}, num2str(jd1));

t0 = cspice_str2et(newstr);
datestring = cspice_et2utc(t0,'C', 6);