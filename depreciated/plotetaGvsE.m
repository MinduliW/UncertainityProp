function plotetaGvsE(nE, nG,param)

labels = {' \eta (\nu)', '\eta(p_1) ', '\eta(p_2) ', '\eta(L) ', '\eta(q_1) ', '\eta(q_2) '};
figure; 
for i = 1:6
    subplot(2,3,i); hold on; 
    plot( param.tvec(1:end-1)*param.TU/60/60, nG(:,i));
    p = plot( param.tvec(1:end-1)*param.TU/60/60, nE(:,i),'--');
    plot_latex(p, 'Time(h)', labels{i},'', '' ,{'Geqoe', 'Eqoe'})

end


end