% check the platform cause my code is now cross platform/annoying 
if ismac
    param.mac = true;
elseif ispc
    param.mac = false;
end

if param.mac == true
    addpath('/Users/minduli/mosek/10.0/toolbox/r2017a');
    addpath('/Users/minduli/libraries_mice/mice/src/mice/');
    addpath('/Users/minduli/libraries_mice/mice/lib');
    addpath('SGP4routines_NAIF')
else
    addpath('H:\Libraries\Mosek\10.0\toolbox\r2017a');
    addpath('H:/Libraries/libraries_mice/mice/src/mice/');
    addpath('H:/Libraries/libraries_mice/mice/lib')
end