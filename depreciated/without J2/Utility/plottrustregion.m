function plottrustregion(param)

fid = fopen("dxmax.txt");
tline = fgetl(fid);

M = [];
while ischar(tline)
    
    % get the second quantity
    newStr = str2double(split(tline));
    arr = newStr(1:end-1)';
    
    M = [M; arr];
    tline = fgetl(fid);
   
   
    
end

[l,c] = size(M);

title = strcat('Trend for nu = ', num2str(param.nu));

figure;  subplot(2,1,1);hold on;
plot(1:l, M(:,1)*param.LU/1e3);
plot(1:l, M(:,2)*param.LU/1e3);
p =plot(1:l, M(:,3)*param.LU/1e3);
plot_latex(p, 'node', 'dx (km)','', title ,{'$\Delta x$', '$\Delta y$', '$\Delta z$'})

subplot(2,1,2); hold on;
plot(1:l, M(:,4)*param.LU/1e3/param.TU);
plot(1:l, M(:,5)*param.LU/1e3/param.TU);
p =plot(1:l, M(:,6)*param.LU/1e3/param.TU);
plot_latex(p, 'node', 'dv (km/s)','', '' ,{'$\Delta v_x$','$\Delta v_y$','$\Delta v_z$'})
end
