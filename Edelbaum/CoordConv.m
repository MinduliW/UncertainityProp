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
            %  coe(4) = right ascension of ascending node (radians)
            %  coe(5) = argument of periapsis (radians)
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

    
    
end
end
