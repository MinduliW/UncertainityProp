function MAT = ECI2RTN(rr,vv)

rhat = rr/norm(rr);
Nhat = cross(rr,vv)/norm(cross(rr,vv));
that = cross(Nhat,rhat);

MAT = [rhat'; that'; Nhat'];
end
