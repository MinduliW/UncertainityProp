% Function turns the figure background white, increases font sizes, and
% allows mathematical formulae be inputted as figure titles. 
% (Makes figures appropriate for a Latex report)
function plot_latex(plt, xtitle, ytitle,ztitle, tl ,str)
% set the plot legend 
if isempty(str) == false 
    hLegend = legend(str, 'Interpreter','latex','Location','NorthOutside','Orientation', 'Horizontal'); 
    hLegend.FontSize = 18;
end 

% set x label, y label and title of the graph, as interpreted using latex 
hXLabel = xlabel(xtitle,'fontsize',20,'Interpreter','latex');
hYLabel= ylabel(ytitle,'fontsize',20,'Interpreter','latex');
zlabel(ztitle,'fontsize',20,'Interpreter','latex');
hTitle = title(tl,'fontsize',20,'Interpreter','latex');

% sets the fontsize of the values on the axes.
set(gca,'fontsize',20)

% the default matlab plot background is grey, this turns it to white 
set(gcf,'color','w');


% increase the width of lines plotted on the graph
% set(plt , 'linewidth',2);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

% Adjust font
set(gca, 'FontName', 'Helvetica')
set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')

% Adjust axes properties
% set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
%     'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
%     'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
%     'LineWidth', 1)

grid minor;
% fullscreen 
% figure('units','normalized','outerposition',[0 0 1 1])

end 
