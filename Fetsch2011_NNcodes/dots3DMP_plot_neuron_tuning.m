function dots3DMP_plot_neuron_tuning(meanFRs,semFRs,u,cohs,hdgs)

tempFR  = meanFRs(:,:,:,:,u);
tempsem = semFRs(:,:,:,:,u);
    
figure('position',[100 500 600 150],'color','w')

deg = char(176);

ax=subplot(1,length(cohs)+1,1); hold on; % single cues
errorbar(hdgs,squeeze(tempFR(1,1,2,:)),squeeze(tempsem(1,1,2,:)),'color','k','marker','o','markerfacecolor','k')
errorbar(hdgs,squeeze(tempFR(2,1,2,:)),squeeze(tempsem(2,1,2,:)),'color','m','marker','v','markerfacecolor','m')
errorbar(hdgs,squeeze(tempFR(2,2,2,:)),squeeze(tempsem(2,2,2,:)),'color','r','marker','s','markerfacecolor','r')
niceformat(ax);
ax.YLim = [0 max(tempFR(:))*1.01];
lh = legend('Vestibular', ['Visual: ' num2str(cohs(1)) '%'], ['Visual: ' num2str(cohs(2)) '%']);
lh.Box = 'off';
title('Single-cue conditions')
ylabel('Firing rate (spikes per s)')

for c=1:length(cohs)
    ax=subplot(1,length(cohs)+1,1+c); hold on
    errorbar(hdgs,squeeze(tempFR(3,c,1,:)),squeeze(tempsem(3,c,1,:)),'color','b','marker','o','markerfacecolor','b')
    errorbar(hdgs,squeeze(tempFR(3,c,2,:)),squeeze(tempsem(3,c,2,:)),'color',rgb('DarkOrange'),'marker','o','markerfacecolor',rgb('DarkOrange'))
    errorbar(hdgs,squeeze(tempFR(3,c,3,:)),squeeze(tempsem(3,c,3,:)),'color',rgb('Green'),'marker','o','markerfacecolor',rgb('Green'))
    niceformat(ax);
    ax.YLim = [0 max(tempFR(:))*1.01];
    ax.YTickLabels = '';
    title(sprintf('Combined, %.d%% coh',cohs(c)))
    if c==1
        xlabel(sprintf('Heading (%s)',deg)); 
        lh=legend(['\Delta = -4' deg],['\Delta = 0' deg],['\Delta = +4' deg]);
        lh.Box = 'off';
    end
end


function niceformat(ax)
ax.TickDir = 'out';
ax.XTick = -10:5:10;
ax.XMinorTick = 'on';
ax.XAxis.MinorTickValues = -10:2.5:10;
ax.YTick = 0:30:100;
ax.YMinorTick = 'on';
ax.YAxis.MinorTickValues = 0:15:100;
ax.TickLength = [0.05 0.05];
ax.XLim = [-11 11];
