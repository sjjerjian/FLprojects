% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

         %ves %vis %comb
clr{1} = {[0 0 0], [255 0 255]./255, [0 255 255]./255};
clr{2} = {[0 0 0], [255 0 0]./255, [0 85 255]./255};
clr{3} = {[0 0 0], [155 0 0]./255, [0 0 225]./255};

modtxt = {'ves','vis','comb'};

for c = 1:length(cohs)
    clear L h;
    figure(c); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    for m = 1:length(mods)     % m c d h
        if m==1
            I = data.modality==mods(m);
        else
            I = data.modality==mods(m) & data.coherence==cohs(c);
        end
        [beta,fval,exitflag,output] = fminsearch(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), [0 3]);
        h(m) = plot(xVals, cgauss(beta,xVals), '-', 'Color', clr{c}{m}, 'Linewidth', 3); hold on;
        errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), 'o', 'Color', clr{c}{m}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        L{m} = modtxt{m};
        xlabel('Heading angle (deg)'); ylabel('Proportion rightward choices'); ylim([0 1]);
        changeAxesFontSize(gca,20,20); set(gca,'box','off')
        legend(h,L,'location','northwest'); legend('boxoff');
        export_fig(['c=' num2str(cohs(c)) '_' modtxt{m} '_pmf'],'-eps');
    end

    figure(100+c); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    for m = 1:length(mods)        
        if m==1
            I = data.modality==mods(m);
        else
            I = data.modality==mods(m) & data.coherence==cohs(c);
        end
        [beta,fval,exitflag,output] = fminsearch(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), [0.7 0 4 0.1]);
        h(d) = plot(xVals, flippedGauss(beta,xVals), '-', 'Color', clr{c}{m}, 'Linewidth', 3); hold on;
        errorbar(hdgs, squeeze(confMean(m,c,D,:)), squeeze(confSE(m,c,D,:)), 'o', 'Color', clr{c}{m}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        xlabel('Heading angle (deg)'); ylabel('Confidence'); ylim([0 1]);
        changeAxesFontSize(gca,20,20); set(gca,'box','off')
        export_fig(['c=' num2str(cohs(c)) '_' modtxt{m} '_conf'],'-eps');
    end

end


% if sum(~isnan(data.rt))>10
%     figure(600+c); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
%     for m = 1:length(mods)        
%         if m==1
%             I = data.modality==mods(m);
%         else
%             I = data.modality==mods(m) & data.coherence==cohs(c);
%         end
%         %         [beta,fval,exitflag,output] = fminsearch(@(x) flippedGauss_err(x,data.rt(I),data.heading(I)), [0.7 0 4 0.1]);
% %         h(d) = plot(xVals, flippedGauss(beta,xVals), '-', 'Color', clr{c}{m}, 'Linewidth', 3); hold on;
%         errorbar(hdgs, squeeze(RTmean(m,c,D,:)), squeeze(RTse(m,c,D,:)), 'o-', 'Color', clr{c}{m}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
%         xlabel('Heading angle (deg)'); ylabel('Reaction time (s)'); ylim([1 2.5]);
%         changeAxesFontSize(gca,20,20); set(gca,'box','off')
%         export_fig(['c=' num2str(cohs(c)) '_' modtxt{m} '_RT'],'-eps');
%     end
% end


%% now separate by delta

if length(deltas)>1


close all

clr{1} = {[255 140 0]./255, [175 238 238]./255, [221 160 221]./255};
clr{2} = clr{1};
clr{3} = clr{1};

for c = 1:length(cohs)
    clear L h;
    figure(200+c); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    for d = 1:length(deltas)     % m c d h
        I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);
        [beta,fval,exitflag,output] = fminsearch(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), [0 3]);
% % %         muPMF(3,c,d) = beta(1);
% % %         sigmaPMF(3,c,d) = beta(2);        
        h(d) = plot(xVals, cgauss(beta,xVals), '-', 'Color', clr{c}{d}, 'Linewidth', 3); hold on;
        errorbar(hdgs, squeeze(pRight(3,c,d,:)), squeeze(pRightSE(3,c,d,:)), 'o', 'Color', clr{c}{d}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        L{d} = sprintf('\\Delta=%d',deltas(d));
        xlabel('Heading angle (deg)'); ylabel('Proportion rightward choices'); ylim([0 1]);
        changeAxesFontSize(gca,20,20); set(gca,'box','off')
        legend(h,L,'location','northwest'); legend('boxoff');
        export_fig(['c=' num2str(cohs(c)) '_d=' num2str(deltas(d)) '_pmf'],'-eps');
    end

    figure(300+c); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    for d = 1:length(deltas)
        I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);
        [beta,fval,exitflag,output] = fminsearch(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), [0.7 0 4 0.1]);
% % %         amplConf(3,c,d) = beta(1);
% % %         muConf(3,c,d) = beta(2);
% % %         sigmaConf(3,c,d) = beta(3);
        h(d) = plot(xVals, flippedGauss(beta,xVals), '-', 'Color', clr{c}{d}, 'Linewidth', 3); hold on;
        errorbar(hdgs, squeeze(confMean(3,c,d,:)), squeeze(confSE(3,c,d,:)), 'o', 'Color', clr{c}{d}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        xlabel('Heading angle (deg)'); ylabel('Confidence'); ylim([0 1]);
        changeAxesFontSize(gca,20,20); set(gca,'box','off')
        export_fig(['c=' num2str(cohs(c)) '_d=' num2str(deltas(d)) '_conf'],'-eps');
    end
end

end