function fx = reconstructPolynomial(filename, dx)


[DAx] = LoadCOSY(filename,6,6,0);

for i = 1: length(DAx)

    DAx(i).C

    constantpart =  DAx(i).C(1); 
    
    for j = 2:length(DAx(i).C(1))

        if sum(DAx(i).E(j,:)) == 1
         additions = additions + dx*DAx(i).E(j,:)*DAx(i).C(j);
        end


    end


    fx(i) = constantpart + additions;



end

end
