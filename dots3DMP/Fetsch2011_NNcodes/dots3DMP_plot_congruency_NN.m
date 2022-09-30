function [congruencyIndex, isSignificant] = dots3DMP_plot_congruency_NN(...
    monkUnit,tuning_curves,vesCorr,visCorr,vesP,visP,xq)
      
% use 60% (high coh), as in paper
congruencyIndex = vesCorr .* visCorr(:,2);
isSignificant = (vesP<0.05) & (visP(:,2)<0.05);


bw = 0.1; % binwidth

figure('color','w','position',[600 200 300 450])
ax=subplot(311); hold on; box off;
histogram(congruencyIndex(monkUnit=='Y'),'BinWidth',bw,'facecolor','w')
histogram(congruencyIndex(monkUnit=='Y' & isSignificant),'BinWidth',bw,'facecolor','k')
axis([-1 1 0 10]);
ax.TickDir = 'out'; 
ax.XTick = -1:0.25:1; ax.XTickLabel = {'-1','','-0.5','','0','','0.5','','1'};
ax.YTick = 0:2.5:10; ax.YTickLabel = {'0','','5','','10'};
title(sprintf('Monkey Y (N = %d)',sum(monkUnit=='Y')));

ax=subplot(312); hold on; box off;
histogram(congruencyIndex(monkUnit=='W'),'BinWidth',bw,'facecolor','w')
histogram(congruencyIndex(monkUnit=='W' & isSignificant),'BinWidth',bw,'facecolor','k')
axis([-1 1 0 6])
ax.TickDir = 'out'; 
ax.XTick = -1:0.25:1; ax.XTickLabel = {'-1','','-0.5','','0','','0.5','','1'};
ax.YTick = 0:1.5:6; ax.YTickLabel = {'0','','3','','6'};
ylabel('Number of neurons')
title(sprintf('Monkey W (N = %d)',sum(monkUnit=='W')));

ax=subplot(313); hold on; box off;
histogram(congruencyIndex,'BinWidth',bw,'facecolor','w')
histogram(congruencyIndex(isSignificant),'BinWidth',bw,'facecolor','k')
axis([-1 1 0 16])
ax.TickDir = 'out'; 
ax.XTick = -1:0.25:1; ax.XTickLabel = {'-1','','-0.5','','0','','0.5','','1'};
ax.YTick = 0:4:16; ax.YTickLabel = {'0','','8','','16'};
title(sprintf('Both monkeys (N = %d)',length(congruencyIndex)));

%%
% flip some tuning curves so that all cells prefer rightward at high vis coh!
isVisRight = visCorr(:,2)>0;
tuning_curves(:,:,:,:,~isVisRight) = tuning_curves(:,:,:,end:-1:1,~isVisRight);

isCongruent = congruencyIndex > 0.4;
isOpposite  = congruencyIndex < -0.4;

figure('color','w','position',[500 500 250 350])
subplot(211); hold on;
plot(xq,squeeze(mean(tuning_curves(1,1,2,:,isCongruent),5)),'color','k','linew',1.5);
plot(xq,squeeze(mean(tuning_curves(2,1,2,:,isCongruent),5)),'color','m','linew',1.5,'linestyle','--');
plot(xq,squeeze(mean(tuning_curves(3,1,2,:,isCongruent),5)),'color','c','linew',1.5,'linestyle','--');
plot(xq,squeeze(mean(tuning_curves(2,2,2,:,isCongruent),5)),'color','r','linew',1.5);
plot(xq,squeeze(mean(tuning_curves(3,2,2,:,isCongruent),5)),'color','b','linew',1.5);
axis([-11 11 10 40])
title('Congruent cells (CI > 0.4)')

subplot(212); hold on;
plot(xq,squeeze(mean(tuning_curves(1,1,2,:,isOpposite),5)),'color','k','linew',1.5);
plot(xq,squeeze(mean(tuning_curves(2,1,2,:,isOpposite),5)),'color','m','linew',1.5,'linestyle','--');
plot(xq,squeeze(mean(tuning_curves(3,1,2,:,isOpposite),5)),'color','c','linew',1.5,'linestyle','--');
plot(xq,squeeze(mean(tuning_curves(2,2,2,:,isOpposite),5)),'color','r','linew',1.5);
plot(xq,squeeze(mean(tuning_curves(3,2,2,:,isOpposite),5)),'color','b','linew',1.5);
axis([-11 11 10 40])
title('Opposite cells (CI < -0.4)')
xlabel(sprintf('Heading (%s)',char(176))); 