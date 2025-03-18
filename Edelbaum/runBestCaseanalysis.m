function runBestCaseanalysis(classes, classnumbers, debrisno, param)

figure;
if debrisno == 3
    t = tiledlayout(3,3);
else
    t= tiledlayout(3,3);
end

if param.Topt == false
    t= tiledlayout(3,3);
end

t.Title.String = strcat('DR = ', {' '}, num2str(param.dutyRatio), {' '}, ', Tmax = ', ...
    {' '}, num2str(param.T), 'N', {' '}, ', Isp = ', num2str(param.Isp), 's');
t.Title.FontWeight = 'bold';
t.Title.FontSize = 20;

for classno = classnumbers
    
    as = classes(classno).as;
    incs = classes(classno).incs;
    
    medianinc = mean(incs);
    deltainc = linspace(0, 0.5*pi/180, 10);
    
    deltaa = linspace(0, 50e3, 10);
    
   
    if param.consta == true
        if debrisno == 3
            [as, deltainc] = meshgrid(as(1:end), deltainc(1:end));
            
        else
            [as, deltainc] = meshgrid(as(1:end), deltainc(1:end));
        end
    else
        
        if debrisno == 3
            [as, deltaa] = meshgrid(as(1:end), deltaa(1:end));
            
        else
            [as, deltaa] = meshgrid(as(1:end), deltaa(1:end));
        end
        
    end
    
    
    
    [rows,cols] = size(as);
    
    
    
    smac = [];incc = [];dvs = []; omegaopt1 = [];
    omegaopt2 = [];mp= [];
    
    
    
    for i = 1:rows
        for j = 1:cols
            
            
            if debrisno == 3
                m_debris = [3e3, 3e3, 3e3];
                
                if param.consta == true
                    debris = [as(i,j), medianinc, 0;...
                        as(i,j), medianinc+deltainc(i,j),0;...
                        as(i,j), medianinc-deltainc(i,j),0];
                else
                    debris = [as(i,j), medianinc, 0;...
                        as(i,j)+deltaa(i,j), medianinc,0;...
                        as(i,j)-deltaa(i,j), medianinc,0];
                end
                
                perm = perms([1,2,3]);
                
            elseif debrisno == 2
                if param.consta == true
                    debris = [as(i,j), medianinc,0;as(i,j), medianinc+2*deltainc(i,j),0];
                else
                    debris = [as(i,j), medianinc,0;as(i,j)+2*deltaa(i,j), medianinc,0];
                    
                end
                
                m_debris = [3e3, 3e3];
                perm = perms([1,2]);
                
            end
            
            
            param.plots = false;
            
            fval = 1e10;
            ddata = [];
            dorder = [];
            for l = 1 %:factorial(debrisno)
                if debrisno == 3
                    debrisorder = [debris(perm(l,1), : );...
                        debris(perm(l,2), : ); debris(perm(l,3), : )];
                else
                    debrisorder = [debris(perm(l,1), : );...
                        debris(perm(l,2), : )];
                end
                
                [RAANs,fvalnew,dDatanew,propmassnew]  = estimOptRAAN(debrisorder,m_debris, param);
                
                if propmassnew ==0
                    'here'
                end
                
                if fvalnew < fval
                    fval = fvalnew;
                    ddata = dDatanew;
                    dorder = debrisorder;
                    propmass = propmassnew;
                end
                
                
            end
            
            debrisData = ddata;
            
            param.krange = -15:1:15;
            
            [result,totalSmaChange, totalIncChange] = ...
                tourGeneral(debrisData, dorder,m_debris, param);
            
            while abs(fvalnew -result)>1
                
                param.krange = -50:1:50;
                [result,totalSmaChange, totalIncChange] = ...
                    tourGeneral(debrisData, dorder,m_debris, param);
            end
            
            
            smac(i,j) = totalSmaChange;
            incc(i,j) = totalIncChange;
            omegaopt1(i,j) = (RAANs(1));
            mp(i,j) = propmass;
            if debrisno == 3
                omegaopt2(i,j) = (RAANs(2));
            end
            
            dvs(i,j) = result;
            
        end
        
    end
    
    if classno == 1
        classtitle = 'LALI';
    elseif classno == 2
        classtitle = 'MALI';
    elseif classno == 3
        classtitle = 'HALI';
    elseif classno == 4
        classtitle = 'LAMI';
    elseif classno == 5
        classtitle = 'MAMI';
    elseif classno == 6
        classtitle = 'HAMI';
    elseif classno == 7
        classtitle = 'LAHI';
    elseif classno == 8
        classtitle = 'MAHI';
    elseif classno == 9
        classtitle = 'HAHI';
        
        
    end
    
    
    if param.Topt == true
        
        if param.consta == true
            nexttile;
            hold on;
            p = contourf(as/1000,deltainc*180/pi,  dvs);
            
            hc = colorbar('EastOutside');
            
            title(hc,'TOF (days)');
           % caxis([560,1200]);
            titlest = strcat( classtitle, ': $i_{ref}$ = ', num2str(medianinc*180/pi),{' ' }, 'deg');
            plot_latex(p, '$a$  (km)', '$\delta i$  (deg)','', ...
                titlest ,{})
        else
            nexttile;
            hold on;
            p = contourf(as/1000,deltaa/1e3,  dvs, linspace(800,1600,50));
            
            hc = colorbar('EastOutside');
            
            title(hc,'TOF (days)');
            caxis([800,1600]);
            titlest = strcat( classtitle, ': TOF (days) (',...
                num2str(debrisno), {' '}, 'debris)');
            
            titlest = strcat( classtitle, ': $i_{ref}$ = ', num2str(medianinc*180/pi),{' ' }, 'deg');
            
            plot_latex(p, '$a$  (km)', '$\delta a$  (km)','', ...
                titlest ,{})
            
        end
        
        
        %         nexttile;% subplot(3,1,2);
        %         hold on;
        %         p = contourf(as/1000,incs*180/pi,  ((omegaopt1))*180/pi);
        %         colorbar('EastOutside');
        %         plot_latex(p, '$a$  (km)', '','', ...
        %             ' RAAN1 (deg)' ,{})
        %
        %         if debrisno == 3
        %             nexttile; %subplot(3,1,3);
        %             hold on;
        %             p = contourf(as/1000,incs*180/pi,  ((omegaopt2))*180/pi);
        %             colorbar('EastOutside');
        %             plot_latex(p, '$\Delta a$  (km)', '','', ...
        %                 'Optimal RAAN2 (deg)' ,{})
        %         end
        %
    else
        
        if param.consta == true
            nexttile;
            
            
            p = contourf(as/1000,deltainc*180/pi, mp);
            
            hc = colorbar('EastOutside');
            
            title(hc,'mp (kg)');
            
            titlest = strcat( classtitle, ': $i_{ref}$ = ', num2str(medianinc*180/pi),{' ' },  'deg');
            
            plot_latex(p, '$a$  (km)', '$\delta i$  (deg)','', ...
                titlest ,{})
            
            
        else
            
            nexttile;
            hold on;
            p = contourf(as/1000,deltaa/1e3,  dvs*1e3 ,linspace(650,1450,50));
            
            hc = colorbar('EastOutside');
            title(hc,'dv (m/s)');
            caxis([650,1450]);
            
            
            
            
            titlest = strcat( classtitle, ': $i_{ref}$ = ', num2str(medianinc*180/pi),{' ' },'deg');
            
            plot_latex(p, '$a$  (km)', '$\delta a$  (km)','', ...
                titlest ,{})
            
            
            
        end
        
    end
    
end


end
