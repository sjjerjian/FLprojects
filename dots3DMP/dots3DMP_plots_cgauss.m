% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

muPMF = nan(length(mods),length(cohs),length(deltas));
muPMFse = nan(length(mods),length(cohs),length(deltas));
sigmaPMF = nan(length(mods),length(cohs),length(deltas));
sigmaPMFse = nan(length(mods),length(cohs),length(deltas));
muConf = nan(length(mods),length(cohs),length(deltas));
muConfse = nan(length(mods),length(cohs),length(deltas));
sigmaConf = nan(length(mods),length(cohs),length(deltas));
sigmaConfse = nan(length(mods),length(cohs),length(deltas));

         %ves %vis %comb
clr{1} = {'ko','mo','co'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
figure(101+D);
set(gcf,'Color',[1 1 1],'Position',[300 500 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(2,length(cohs),c);
    for m = 1:length(mods)     % m c d h
        if m==1
            I = data.modality==mods(m);
        else
            if D==length(deltas)+1
                I = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
            else
                I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(D);
            end
        end
        beta = fminsearch(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), [0 3]);
        [betaUnc,~,~,~,~,hessian] = fminunc(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), [0 3]);
        SE = sqrt(diag(inv(hessian)));
        muPMFse(m,c,D) = SE(1);
        sigmaPMFse(m,c,D) = SE(2);
        if unc
            muPMF(m,c,D) = betaUnc(1);
            sigmaPMF(m,c,D) = betaUnc(2);
        else
            muPMF(m,c,D) = beta(1);
            sigmaPMF(m,c,D) = beta(2);
        end

        h(m) = plot(xVals, cgauss(beta,xVals), [clr{c}{m}(1) '-']); hold on;
        errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), clr{c}{m});
        ylim([0 1]); if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,'vestib','visual','comb','Location','northwest');
    xlabel('heading angle (deg)'); ylabel('proportion rightward choices');

    subplot(2,length(cohs),c+2+(length(cohs)-2));
    for m = 1:length(mods)        
        if m==1
            I = data.modality==mods(m);
        else
            if D==length(deltas)+1
                I = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
            else
                I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(D);
            end
        end
        beta = fminsearch(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), [0.7 0 4 0.1]);
        [betaUnc,~,~,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), [0.7 0 4 0.1]);
        SE = sqrt(diag(inv(hessian)));
        amplConfse(m,c,D) = SE(1);
        muConfse(m,c,D) = SE(2);
        sigmaConfse(m,c,D) = SE(3);
        if unc
            amplConf(m,c,D) = betaUnc(1);
            muConf(m,c,D) = betaUnc(2);
            sigmaConf(m,c,D) = betaUnc(3);
        else
            amplConf(m,c,D) = beta(1);
            muConf(m,c,D) = beta(2);
            sigmaConf(m,c,D) = beta(3);
        end
        
        h(m) = plot(xVals, flippedGauss(beta,xVals), [clr{c}{m}(1) '-']); hold on;       
        errorbar(hdgs, squeeze(confMean(m,c,D,:)), squeeze(confSE(m,c,D,:)), clr{c}{m});
        ylim([0 1]); if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
%     legend(h,'vestib','visual','comb','location','northwest');
    xlabel('heading angle (deg)'); ylabel('saccadic endpoint (''confidence'', %)');

end


%% now separate by delta

clr{1} = {'bs','cs','gs'};
clr{2} = {'b^','c^','g^'};
clr{3} = {'bo','co','go'};

clear L;
figure(108);
set(gcf,'Color',[1 1 1],'Position',[50 20 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(2,length(cohs),c);
    for d = 1:length(deltas)     % m c d h
        I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);
        beta = fminsearch(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), [0 3]);
        [betaUnc,~,~,~,~,hessian] = fminunc(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), [0 3]);
        SE = sqrt(diag(inv(hessian)));
        muPMFse(3,c,d) = SE(1);
        sigmaPMFse(3,c,d) = SE(2);        
        if unc
            muPMF(3,c,d) = betaUnc(1);
            sigmaPMF(3,c,d) = betaUnc(2);
        else
            muPMF(3,c,d) = beta(1);
            sigmaPMF(3,c,d) = beta(2);
        end
        
        h(d) = plot(xVals, cgauss(beta,xVals), [clr{c}{d}(1) '-']); hold on;
        errorbar(hdgs, squeeze(pRight(3,c,d,:)), squeeze(pRightSE(3,c,d,:)), clr{c}{d});
        L{d} = sprintf('\\Delta=%d',deltas(d));
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,L,'location','northwest');
    xlabel('heading angle (deg)'); ylabel('proportion rightward choices');

    subplot(2,length(cohs),c+2+(length(cohs)-2));
    for d = 1:length(deltas)
        I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);
        beta = fminsearch(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), [0.7 0 4 0.1]);
        [betaUnc,~,~,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), [0.7 0 4 0.1]);
        SE = sqrt(diag(inv(hessian)));
        amplConfse(3,c,d) = SE(1);
        muConfse(3,c,d) = SE(2);
        sigmaConfse(3,c,d) = SE(3);
        if unc
            amplConf(3,c,d) = betaUnc(1);
            muConf(3,c,d) = betaUnc(2);
            sigmaConf(3,c,d) = betaUnc(3);
        else
            amplConf(3,c,d) = beta(1);
            muConf(3,c,d) = beta(2);
            sigmaConf(3,c,d) = beta(3);
        end
        h(d) = plot(xVals, flippedGauss(beta,xVals), [clr{c}{d}(1) '-']); hold on;       
        errorbar(hdgs, squeeze(confMean(3,c,d,:)), squeeze(confSE(3,c,d,:)), clr{c}{d});
        L{d} = sprintf('?=%d',deltas(d));
        ylim([0 1]); hold on;
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
%     legend(h,L,'location','northwest');
    xlabel('heading angle (deg)'); ylabel('saccadic endpoint (''confidence'', %)');
end