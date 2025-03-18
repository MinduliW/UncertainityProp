function solveWithDrag2(waitTime,ad,Id, ...
    RAAN0,massAfterL1, timeAfterL1, param)

% output : residual on final position acheived.
% input : wait time

accel_drag = drag_acceleration(ad, timeAfterL1+waitTime, param);


end
