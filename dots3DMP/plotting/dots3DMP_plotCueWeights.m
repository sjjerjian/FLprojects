function dots3DMP_plotCueWeights(wves,wvesBoot,cohs,conftask)

figure(810); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
er = errorbar(cohs,wves.choice.pred,std(wvesBoot.choice.pred),'k-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','k'); hold on;
set(er,'Color','k');
ylim([0 1.1]);
set(gca,'XTick',cohs,'Xlim',[cohs(1)-0.04 cohs(end)+0.04],'XTickLabel',cohs,'Ytick',0:0.2:1);
xlabel('Visual coherence'); ylabel('Vestibular weight');
changeAxesFontSize(gca,20,20); set(gca,'box','off');
% l = legend('Predicted (from single-cues)'); legend('boxoff');
% set(l,'Position',[0.2484    0.9250    0.7444    0.0812]);

er = errorbar(cohs,wves.choice.emp,std(wvesBoot.choice.emp),'m-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','m');    
set(er,'Color','m')
% l = legend('Predicted', 'Empirical');
% set(l,'Position',[0.2317    0.8492    0.7444    0.1516]);

if conftask
    er = errorbar(cohs,wves.conf.emp,std(wvesBoot.conf.emp),'c-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','c');
    set(er,'Color','c')
    l = legend('Predicted', 'Empirical', 'From confidence curves'); legend('boxoff');
    set(l,'Position',[0.3095    0.7890    0.7444    0.2219]);
end
