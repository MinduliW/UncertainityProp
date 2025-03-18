function deltav = getDvDirect(x, param)


Vd = x(1); %m/s
Id = x(2); % rad


a0 = param.mu/param.x0(1)^2; %m
ad = param.mu/Vd^2; %m
af = param.mu/param.xf(1)^2; %m

[~,~, ~, ~, deltaV_Edel1,~,~] = Kechichan_Algorithm(a0, ad, param.x0(2),Id,param);
deltaV1 = 0.9766*deltaV_Edel1(end);

[~,~, ~, ~, deltaV_Edel2,~,~] = Kechichan_Algorithm(ad, af, Id ,param.xf(2),param);
deltaV2 = 0.9766*deltaV_Edel2(end);

deltav = deltaV1 + deltaV2;



end
