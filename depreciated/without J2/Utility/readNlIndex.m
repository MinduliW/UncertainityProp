% function to extract the limits of each variable given the nonlinearity
% index.
%
function Morg = readNlIndex(fn)

fid = fopen(fn);
tline = fgetl(fid);
lc = 13;
enc = 0;
CA =[];M = [];
while ischar(tline)
    
    % get the second quantity
    newStr = split(tline);
    
    if (rem(lc,13) ==0)
        
        if ~isempty(CA)
            M = [M;CA];
        end
        
        enc = 1;
        CA =[];
    end
    
    if (enc== 12 || enc == 13 ||  enc == 1)
    else
        CA =  [CA,  str2double(string(newStr(3)))];
    end

    tline = fgetl(fid);
    enc = enc+1;
    
    
    lc = lc+1;
    
    
end
   M = [M;CA];
   
% must reorgnaise this so that it's [xdynamics and then control!]

Morg = [M(:,1), M(:, 6:10), M(:,2:5)];

