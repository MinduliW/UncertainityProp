function write2debris(classno,TOFs,fuelcs,thrusttime)

fid = fopen('2debris.txt','w');
fprintf(fid, '%s\n',num2str(classno),num2str(TOFs),num2str(fuelcs),num2str(thrusttime));
fclose(fid);

end
