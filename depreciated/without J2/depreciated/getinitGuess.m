function [xguess,dvSumprev] =getinitGuess(xoCartOsc, leg1, dVendPrev, param)

xguess = [];
dV = zeros(size(param.tvec));

rs = xoCartOsc;

coe(1,:) = CoordConv.vec2orbElem(rs(1:3),rs(4:6),param.mu);

for j = 1: length(param.tvec)-1

    arglat = wrapTo2Pi(coe(j,5)+coe(j,6));

    % calculate components of dv
    beta = interp1(leg1.t, leg1.beta, param.tvec(j));
    dV(j) = interp1(leg1.t, leg1.dV, param.tvec(j));
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

        dvRTN = (dV(j)-dVendPrev)*[0, cos(beta), sin(beta)];
    end


    % convert dv to ECI frame

     dvECI = ECI2RTN(rs(1:3),rs(4:6))'*dvRTN';
       

     % add dv to the current state
    rs(4:6) = rs(4:6)+dvECI;

 
    % propagate from current node to the next.
    [tProp,propstate] = ode45(@(t,x) propagateCart(t, x, param),...
        [param.tvec(j),param.tvec(j+1)],rs, odeset('RelTol', 1e-13,'AbsTol',1e-13));


    res = [rs(1:3)',rs(4:6)',dvECI'];

    % convert to osculating 
    %oscCart = mean2oscCart(rs(1:3)',rs(4:6)',param.mu, param.J2, param.Re);

    %resOsc = [oscCart',dvECI'] ;
    xguess = [xguess; res];

    % make rs the propagated state.
    rs =propstate(end,:)';

    coe(j+1,:) = CoordConv.vec2orbElem(rs(1:3),rs(4:6),param.mu);


end

%oscCart = mean2oscCart(rs(1:3)',rs(4:6)',param.mu, param.J2, param.Re);

resOsc = [rs',0,0,0] ;
xguess = [xguess; resOsc];

dvSumprev = dV(end-1);


end

