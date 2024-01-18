function dots3DMP_plots_NN(parsedData,mods,cohs,deltas,hdgs)


%% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

         %ves %vis %comb
clr{1} = {'ko','mo','co'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
figure(101+D);
set(gcf,'Color',[1 1 1],'Position',[300 1000 450+300*(length(cohs)-2) 500],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(2,length(cohs),c); box off; hold on;
    for m = 1:length(mods)     % m c d h
        if parsedData.plotLogistic(m,c,D)
            h(m) = plot(parsedData.xVals,squeeze(parsedData.yVals(m,c,D,:)),[clr{c}{m}(1) '-'],'linewidth',1.5); hold on;
            errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:)), squeeze(parsedData.pRightSE(m,c,D,:)), clr{c}{m},'linewidth',1.5);
        else
            h(m) = errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:)), squeeze(parsedData.pRightSE(m,c,D,:)), [clr{c}{m} '-'],'linewidth',1.5); hold on;
        end
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,'vestib','visual','comb','Location','northwest');
    xlabel('heading angle (deg)'); ylabel('proportion rightward choices');
    try, changeAxesFontSize(gca,15,15); end
end

set(gcf,'Color',[1 1 1],'Position',[300 1000 450+300*(length(cohs)-2) 500],'PaperPositionMode','auto'); clf;


%% now separate by delta

if length(deltas)>1

clr{1} = {'bs','cs','gs'};
clr{2} = {'b^','c^','g^'};
clr{3} = {'bo','co','go'};

clear L;
%figure(108);
%set(gcf,'Color',[1 1 1],'Position',[50 20 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(2,length(cohs),c+length(cohs)); box off; hold on;
    for d = 1:length(deltas)     % m c d h
        if parsedData.plotLogistic(3,c,d)
            h(d) = plot(parsedData.xVals,squeeze(parsedData.yVals(3,c,d,:)),[clr{c}{d}(1) '-'], 'linewidth', 1.5); hold on;
            errorbar(hdgs, squeeze(parsedData.pRight(3,c,d,:)), squeeze(parsedData.pRightSE(3,c,d,:)), clr{c}{d});
        else
            h(d) = errorbar(hdgs, squeeze(parsedData.pRight(3,c,d,:)), squeeze(parsedData.pRightSE(3,c,d,:)), [clr{c}{d} '-'], 'linewidth', 1.5); hold on
        end
        L{d} = sprintf('\x0394=%d',deltas(d));
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,L,'location','northwest');
    xlabel('heading angle (deg)'); ylabel('proportion rightward choices');
    try, changeAxesFontSize(gca,15,15); end
end

end