function [c,ceq] = contraints(x, debris,m_debris, param)

[~,~, ~,Total_TOF, Total_dV,Total_fuel,~,decaystatus] = ...
    tourGeneral(x, debris,m_debris, param);




if param.Topt == true

  %c= Total_dV - param.dvLimit; 
  c(1) = Total_fuel - param.fuellimit; %kg
  %decaystates <1 
  %c(2) = decaystatus -1; 
else
%    result = stepfcn(max(wt), param.waitTimeLimit)*(max(wt) -...
%         param.waitTimeLimit)^2 + Total_dV;  
  c(1) = Total_TOF - param.toflimit*param.TU; 
 % c(2) = decaystatus -1; 
end

ceq = []; 


end
