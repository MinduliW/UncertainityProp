function  plotTR(dynamics,TR,param)


labels = {'\nu', ' p_1 ', ' p_2 ', ' L ', ' q_1 ', ' q_2 '};
figure; 
for i = 1:6
    subplot(3,2,i); hold on; 
    plot( param.tvec*param.TU/60/60, dynamics(i,:));
    plot( param.tvec(1:end-1)*param.TU/60/60, dynamics(i,(1:end-1))' + TR(:,i),'r');
    p = plot( param.tvec(1:end-1)*param.TU/60/60, dynamics(i,(1:end-1))' - TR(:,i), 'r');
    plot_latex(p, 'Time(h)', labels{i},'', '' ,{})

end

end
