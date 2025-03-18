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
            
            e = abs(e);
            a = abs(a);
            %ecentric anomaly
            cosE = (1 - r./a)./e;
            
            sinE = sum(rs.*vs)./(e.*sqrt(mus.*a));
            E = atan2(sinE,cosE);
            
            %inclination
            sinI = sqrt(L(1,:).^2 + L(2,:).^2)./sqrt(L2s);
            cosI = L(3,:)./sqrt(L2s);
            I = atan2(sinI,cosI);
            
            sinw = ((vs(1,:).*L(2,:) - vs(2,:).*L(1,:))./mus - ...
                rs(3,:)./r)./(e.*sinI);
            cosw = ((sqrt(L2s).*vs(3,:))./mus - (L(1,:).*rs(2,:) - ...
                L(2,:).*rs(1,:))./(sqrt(L2s).*r))./(e.*sinI);
            RAAN = atan2(sinw,cosw);
         
            
            cosO = -L(2,:)./(sqrt(L2s).*sinI);
            sinO = L(1,:)./(sqrt(L2s).*sinI);
            AOP = atan2(sinO,cosO);
            
            if e  <0
                'here'
            end

            %orbital periods
            P = 2*pi*sqrt(a.^3./mus);
            
            %time of periapsis crossing
            tau = -(E - e.*sin(E))./sqrt(mus.*a.^-3);
            
            
            True_an = 2 * atan(sqrt((1 + e)/ (1 - e)) .* tan(E/2) );
            
            if dot(rs,vs) < 0

                True_an = wrapTo2Pi(True_an);

                if True_an > pi && True_an < 2*pi
                else
                    True_an = 2.0*pi - True_an;
                end

               
            end

            
            %time of periapsis crossing
            
            tau = -(E - e.*sin(E))./sqrt(mus.*a.^-3);
            
            arglat = wrapTo2Pi(wrapTo2Pi(True_an) + wrapTo2Pi(AOP)); 
            
            truelon =  wrapTo2Pi(wrapTo2Pi(True_an) + wrapTo2Pi(AOP) + wrapTo2Pi(RAAN));

            x = [a,e,I,AOP, RAAN,True_an,truelon];
            
            
            %A and B vectors
            A = [a.*(cos(AOP).*cos(RAAN) - sin(AOP).*cos(I).*sin(RAAN));...
                a.*(sin(AOP).*cos(RAAN) + cos(AOP).*cos(I).*sin(RAAN));...
                a.*sin(I).*sin(RAAN)];
            B = [-a.*sqrt(1-e.^2).*(cos(AOP).*sin(RAAN) + ...
                sin(AOP).*cos(I).*cos(RAAN));...
                a.*sqrt(1-e.^2).*(-sin(AOP).*sin(RAAN) + ...
                cos(AOP).*cos(I).*cos(RAAN));...
                a.*sqrt(1-e.^2).*sin(I).*cos(RAAN)];
            
            
        end