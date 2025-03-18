function coe = fPropPMDT(tProps,leg1, param)

x0Kep = [leg1.a(1), 0, leg1.inc(1), wrapTo2Pi(leg1.RAAN(1)), 0,wrapTo2Pi(sqrt(param.mu/leg1.a(1)^3)*leg1.t(1))];
% convert osculating starting pos to cartesian
[rr, vv] = CoordConv.po2pv(x0Kep,param.mu);
x0cart =  [rr;vv];

% convert to osculating 


x0CartOsc = mean2oscCart(x0cart(1:3)',x0cart(4:6)',param.mu, param.J2, param.Re);


coe(1,:) = CoordConv.vec2orbElem(x0cart(1:3),x0cart(4:6),param.mu);
coeOsc = mean2oscKep(coe(1,:), param.mu, param.J2, param.Re);

rs = x0CartOsc;


for j = 1: length(tProps)-1

    arglat = wrapTo2Pi(coeOsc(5)+coeOsc(6));

    % calculate components of dv
    beta = interp1(leg1.t, leg1.beta, tProps(j));
    dV(j) = interp1(leg1.t, leg1.dV, tProps(j));
    if  param.incincrease == true
        if (arglat>pi/2 && arglat<=3*pi/2)
            beta = -abs(beta);

        else
            beta = abs(beta);
        end
    else
        if (arglat>=pi/2 && arglat<3*pi/2)
            beta = abs(beta);
        else
            beta = -abs(beta);
        end

    end

    if j > 1
        dvRTN = (dV(j)-dV(j-1))*[0, cos(beta), sin(beta)];
    else

        dvRTN = (dV(j))*[0, cos(beta), sin(beta)];
    end


    % convert dv to ECI frame

     dvECI = ECI2RTN(rs(1:3),rs(4:6))'*dvRTN';
       

     % add dv to the current state
    rs(4:6) = rs(4:6)+dvECI;

    if tProps(j) == tProps(j+1)
        tProp = tProps(j);
        propstate = rs';
    else
        % propagate from current node to the next.
        [tProp,propstate] = ode45(@(t,x) propagateCart(t, x, param),...
            [tProps(j),tProps(j+1)],rs, odeset('RelTol', 1e-13,'AbsTol',1e-13));



    end

    % make rs the propagated state.
    rs =propstate(end,:)';


      % convert to osculating 
      meanCart = osc2meanCart(rs(1:3)',rs(4:6)',param.mu, param.J2, param.Re);
      
      coeOsc = CoordConv.vec2orbElem(propstate(end,1:3),propstate(end,4:6),param.mu);
   
    %rs = oscCart;
    coe(j+1,:) = CoordConv.vec2orbElem(meanCart(1:3),meanCart(4:6),param.mu);


end



end
