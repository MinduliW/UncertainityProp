%Coordinate conversions
classdef CoordConv
    methods (Static)
        function MEEParameters = kepler2MEOE(Object)
            % input: A struct or array with keplerian orbital elements
            % output: A vector of MEE elements.
            
            if isstruct(Object)
                
                a = Object.a;
                e = Object.e;
                omega = Object.omega;
                Omega = Object.Omega;
                i = Object.i;
                
                if isfield(Object, 'MA')
                    MA = Object.MA;
                    
                    
                    E = MA;
                    E_old = 1;
                    precision  = 1e-7;
                    
                    
                    while abs(E - E_old) > precision
                        
                        % find the next eccentric anomaly value using N-R method
                        E = E - ((E - e*sin(E) - MA)./(1 - e*cos(E)));
                        E_old = E ;
                    end
                    
                    theta = 2 * atan(sqrt((1 + e)/ (1 - e)) .* tan(E/2) );
                else
                    theta = Object.theta;
                end
                
                
            else
                a = Object(1);
                e = Object(2);
                omega = Object(3);
                Omega = Object(4);
                i = Object(5);
                theta = Object(6);
                
            end
            
            p = a*(1 - e*e);
            f = e*cos(omega + Omega);
            g = e*sin(omega + Omega);
            h = tan(i/2)*cos(Omega);
            k = tan(i/2)*sin(Omega);
            L = wrapTo2Pi(Omega + omega + theta);
            MEEParameters = [p, f,g ,h,k,L];
            
            
            
        end
        function [rr, vv] = po2pv(PO, mu)
            
            % input: PO vector of classical orbital parameters
            %        mu gravity parameter in units consistent with a
            % output: rr position vecotor in units consistent with mu and a
            %         vv velocity vector in units consistent with mu and a
            
            a = PO(1);
            e = PO(2);
            i = PO(3);
            Om = PO(4);
            om = PO(5);
            theta = PO(6);
            
            A = [cos(om+theta) -sin(om+theta) 0;
                sin(om+theta) cos(om+theta)  0;
                0             0   1];
            
            if i<0
                i = pi+i;
            end
            
            B = [1      0       0;
                0 cos(i)  -sin(i);
                0 sin(i)   cos(i)];
            
            C =  [cos(Om) -sin(Om) 0;
                sin(Om) cos(Om)  0;
                0       0   1];
            
            p = a*(1-e^2);
            
            r = [1/(1+e*cos(theta))*p 0 0]';
            v = sqrt(mu/p)*[e*sin(theta) 1+e*cos(theta) 0]';
            
            rr = C*B*A*r;
            vv = C*B*A*v;
        end
        function OPmat = KeplerStruct(EP)
            OP = CoordConv.ep2op(EP);
            OPmat.aAU = OP(1);
            OPmat.e = OP(2);
            OPmat.IncDeg = rad2deg((OP(3)));
            OPmat.OmegaDeg = rad2deg((OP(4)));
            OPmat.omegaDeg = rad2deg((OP(5)));
            OPmat.TrueAnDeg = rad2deg((OP(6)));
            OPmat.ArgofLat = (OPmat.TrueAnDeg + OPmat.omegaDeg);
            
        end
        
        function OPmat = KeplerStruct2(EP)
            OP = CoordConv.mee2coe(EP);
            OPmat.aAU = OP(1);
            OPmat.e = OP(2);
            OPmat.IncDeg = rad2deg((OP(3)));
            OPmat.OmegaDeg = rad2deg((OP(4)));
            OPmat.omegaDeg = rad2deg((OP(5)));
            OPmat.TrueAnDeg = rad2deg((OP(6)));
            OPmat.ArgofLat = (OPmat.TrueAnDeg + OPmat.omegaDeg);
            
        end
        
        function coe = mee2coe(mee)
            
            % convert modified equinoctial elements to classical orbit elements
            
            % input
            
            %  mee(1) = semiparameter (kilometers)
            %  mee(2) = f equinoctial element
            %  mee(3) = g equinoctial element
            %  mee(4) = h equinoctial element
            %  mee(5) = k equinoctial element
            %  mee(6) = true longitude (radians)
            
            % output
            
            %  coe(1) = semimajor axis (kilometers)
            %  coe(2) = eccentricity
            %  coe(3) = inclination (radians)
            %  coe(4) = argument of periapsis (radians)
            %  coe(5) = right ascension of ascending node (radians)
            %  coe(6) = true anomaly (radians)
            
            % Orbital Mechanics with MATLAB
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % unload modified equinoctial orbital elements
            
            pmee = mee(1);
            fmee = mee(2);
            gmee = mee(3);
            hmee = mee(4);
            kmee = mee(5);
            lmee = mee(6);
            
            % compute classical orbital elements
            
            tani2s = sqrt(hmee * hmee + kmee * kmee);
            
            % orbital eccentricity
            
            ecc = sqrt(fmee * fmee + gmee * gmee);
            
            % semimajor axis
            
            sma = pmee / (1.0 - ecc * ecc);
            
            % orbital inclination
            
            inc = 2.0 * atan(tani2s);
            
            % right ascension of ascending node
            
            raan = atan2(kmee, hmee);
            
            % argument of periapsis
            
            atopo = atan2(gmee, fmee);
            
            argper = mod(atopo - raan, 2.0 * pi);
            
            % true anomaly
            
            tanom = mod(lmee - atopo, 2.0 * pi);
            
            % load classical orbital element array
            
            coe(1) = sma;
            coe(2) = ecc;
            coe(3) = inc;
            coe(4) = raan;
            coe(5) = argper;
            coe(6) = tanom;
        end
        
        function  OP = ep2op(EP)
            % converts equinoctial parameters in orbital parameters
            % EP = (p,f,g,h,k,L)
            % OP = (a,e,i,Om,om,theta)
            
            % Initialize:
            p = EP(1);
            f = EP(2);
            g = EP(3);
            h = EP(4);
            k = EP(5);
            L = (EP(6));
            
            % % Compute:
            % OP(1) = p/(1-f^2-g^2);
            % OP(2) = sqrt(f^2+g^2);
            % OP(3) = 2*atan(sqrt(h^2+k^2));  % Problems finding right quadrant?
            % OP(4) = atan2(g,f)-atan2(k,h);  % Problems finding right quadrant?
            % OP(5) = atan2(k,h);             % Problems finding right quadrant?
            % OP(6) = L- atan2(g,f);          % Problems finding right quadrant?, check by atan(g/f)=O+w.
            
            OP(1) = p/(1-f^2-g^2);
            OP(2) = sqrt(f^2+g^2);
            OP(3) = atan2(2*sqrt(h^2+k^2), 1-h^2-k^2);
            
            if EP(4)==0&&EP(5)==0
                OP(4) = 0;
            else
                OP(4) = atan2(k,h);
            end
            if EP(2)==0&&EP(3)==0
                OP(5) = 0;
            else
                OP(5) = atan2(g*h -f*k,f*h+g*k);
            end
            OP(6) = L - OP(4) - OP(5);
        end
        
        
        function posandvel = ep2pv(EP, mu)
            % Function converts MEE coordinates to positional coordinates
            OP = CoordConv.ep2op(EP);
            
            [rr, vv] = CoordConv.po2pv(OP, mu);
            posandvel = [rr;vv];
        
        end
        
        function x = vec2kepStruct(rs,vs,mus)
            
            OP = CoordConv.vec2orbElem(rs,vs,mus);
            
            x.a = OP(1);
            x.e = OP(2);
            x.IncDeg = rad2deg((OP(3)));
            x.OmegaDeg = rad2deg((OP(5)));
            x.omegaDeg = rad2deg((OP(4)));
            x.TrueAnDeg = rad2deg((OP(6)));
            x.ArgofLat = (x.TrueAnDeg + x.omegaDeg);
            
        end
        
        function mee = vec2mee(rs,vs,mus)
            
            % Get the keplerian elements
            x = CoordConv.vec2orbElem(rs,vs,mus);
            
            a = x(1);
            e = x(2);
            I = x(3); 
            omega = x(4);
            Omega = x(5);
            True_an = x(6);
            
            % convert kepler to mee
            mee =  CoordConv.kepler2MEOE([a,e, omega,Omega, I, True_an]);
        end
        
        
        function x = vec2orbElem(rs,vs,mus)
            
            mus = mus(:).';
            nplanets = numel(rs)/3;
            if mod(nplanets,1) ~= 0 || numel(vs) ~= nplanets*3 ||...
                    (length(mus) ~= nplanets && length(mus) ~= 1)
                error('vec2orbElem:inputError',['rs and vs must contain 3n ',...
                    'elements and mus must have length n or 1 for n bodies']);
            end
            if length(rs) == numel(rs)
                rs = reshape(rs,3,nplanets);
            end
            if length(vs) == numel(vs)
                vs = reshape(vs,3,nplanets);
            end
            v2s = sum(vs.^2);
            r = sqrt(sum(rs.^2)); %orbital separation
            Ws = 0.5*v2s - mus./r;
            a = -mus/2./Ws; %semi-major axis
            L = [rs(2,:).*vs(3,:) - rs(3,:).*vs(2,:);...
                rs(3,:).*vs(1,:) - rs(1,:).*vs(3,:);...
                rs(1,:).*vs(2,:) - rs(2,:).*vs(1,:)]; %angular momentum
            L2s = sum(L.^2);
            p = L2s./mus; %semi-latus rectum
            e = sqrt(1 - p./a); %eccentricity
            %ecentric anomaly
            cosE = (1 - r./a)./e;
            sinE = sum(rs.*vs)./(e.*sqrt(mus.*a));
            E = atan2(sinE,cosE);
            
            %inclination
            sinI = sqrt(L(1,:).^2 + L(2,:).^2)./sqrt(L2s);
            cosI = L(3,:)./sqrt(L2s);
            I = atan2(sinI,cosI);
            %argument of pericenter
            sinw = ((vs(1,:).*L(2,:) - vs(2,:).*L(1,:))./mus - ...
                rs(3,:)./r)./(e.*sinI);
            cosw = ((sqrt(L2s).*vs(3,:))./mus - (L(1,:).*rs(2,:) - ...
                L(2,:).*rs(1,:))./(sqrt(L2s).*r))./(e.*sinI);
            omega = atan2(sinw,cosw);
            %longitude of ascending node
            cosO = -L(2,:)./(sqrt(L2s).*sinI);
            sinO = L(1,:)./(sqrt(L2s).*sinI);
            Omega = atan2(sinO,cosO);
            
            %orbital periods
            P = 2*pi*sqrt(a.^3./mus);
            
            %time of periapsis crossing
            tau = -(E - e.*sin(E))./sqrt(mus.*a.^-3);
            
            True_an = 2 * atan(sqrt((1 + e)/ (1 - e)) .* tan(E/2) );
            
            
            %time of periapsis crossing
            tau = -(E - e.*sin(E))./sqrt(mus.*a.^-3);
            
            x = [a,e,I,Omega,omega,True_an];
            
            
            %A and B vectors
            A = [a.*(cos(Omega).*cos(omega) - sin(Omega).*cos(I).*sin(omega));...
                a.*(sin(Omega).*cos(omega) + cos(Omega).*cos(I).*sin(omega));...
                a.*sin(I).*sin(omega)];
            B = [-a.*sqrt(1-e.^2).*(cos(Omega).*sin(omega) + ...
                sin(Omega).*cos(I).*cos(omega));...
                a.*sqrt(1-e.^2).*(-sin(Omega).*sin(omega) + ...
                cos(Omega).*cos(I).*cos(omega));...
                a.*sqrt(1-e.^2).*sin(I).*cos(omega)];
            
            
        end
    end
    
    
    
    
    
end
