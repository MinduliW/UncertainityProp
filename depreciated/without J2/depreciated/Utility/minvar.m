% calculate the minimum variation allowed to reach a particular
% nonlinearity index for each element

function [limits] = minvar(nu,N)
M = readNlIndex('nonlinIndex.txt');

limits =[];
for i = 1:9
    for j = 1: N
        
        xb(j) = (nu)/abs(M(j, i+1));
        
        % xb(j) = abs(nu-abs(M(j,i)))/abs(M(j, i+1));
        %xub(j) = (-nu)/abs(M(j, i+1));
        
        
    end
    limits = [limits; xb];
    
    %     figure;
    %     plot(xb)
    %min(x)
    %max(xb)
end
end

