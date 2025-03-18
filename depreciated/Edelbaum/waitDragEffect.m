function [Omega_t2, semiMajor,dv, mf] = waitDragEffect(Omega_t1, ad0, Id0,waittime,...
    mass ,param)


% waittime = 24*60*60*1000;
% 
% ad0 = param.Re+ 300e3;
V0 = sqrt(param.mu/ad0);

Omega_t2 = Omega_t1 -param.k*V0^7*cos(Id0)*(waittime);

% if param.drag == true
    
% segment wait time.
waitvec = linspace(0, waittime, 100);

semiMajor(1) = ad0;

f_d_k = -drag_acceleration(ad0, [], param)*ones(size(waitvec));

dv = trapz(waitvec,abs(f_d_k));

% contribution to fuel mass 
endmass = mass/(exp(dv/param.Isp/param.g0));

mf = mass - endmass; 

