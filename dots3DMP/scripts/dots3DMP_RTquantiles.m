function dots3DMP_RTquantiles(data,conftask,correct)

% correct - use correct trials only (1), incorrect only (0), or all (-1)

if nargin < 3, correct = -1; end % use all
  

% loop over 'conditions' indexed as [modality coherence], with the
% familiar manual kluge necessary to avoid invalid combinations:
ucoh = unique(data.coherence);
ucond = [1 ucoh(1); 2 ucoh(1); 2 ucoh(2); 3 ucoh(1); 3 ucoh(2)];
subplotInd = [2 3 4 5 6 1];
titles = {'Ves';'Vis-lo';'Vis-hi';'Comb-lo';'Comb-hi';'All'};

nbins = 5;
uhdg  = unique(abs(data.heading));

mcols = {'Greys','Reds','Reds','Blues','Blues','Purples'};

for c = 1:size(ucond,1)+1 % the extra one is for all conditions pooled
    
    cmap = flipud(cbrewer('seq',mcols{c},length(uhdg)*2));
    cmap = cmap(end-1:-2:1,:);
    
    for h = 1:length(uhdg)
        if c==size(ucond,1)+1
            I = abs(data.heading)==uhdg(h);
        else
            I = abs(data.heading)==uhdg(h) & data.modality==ucond(c,1) & data.coherence==ucond(c,2);%& data.corr==0;
        end
        if correct>=0 % select only correct/incorrect trials (and all heading==0)
            I = I & (data.correct==correct | data.heading==0); 
        else
            temp = data.correct;
            theseCorr = temp(I);
        end
        theseRT = data.RT(I);
        theseConf = data.PDW(I);
        
        rtQ = [0 quantile(theseRT,nbins-1) inf]; % or quintiles
        for q = 1:length(rtQ)-1
            J = theseRT>=rtQ(q) & theseRT<rtQ(q+1);
            X(c,h,q) = mean(theseRT(J));
            Y(c,h,q) = mean(theseConf(J));
            Ye(c,h,q) = std(theseConf(J)) / sqrt(sum(J));
            
            if correct<0
                Yc(c,h,q) = mean(theseCorr(J));
                Yce(c,h,q) = sqrt( (Yc(c,h,q).*(1-Yc(c,h,q))) ./ sum(theseCorr(J)) );
            end
        end
    end

    figure(16+correct);
    set(gcf,'Color',[1 1 1],'Position',[200 200 270*2 170*3],'PaperPositionMode','auto');
    subplot(3,2,subplotInd(c));

    if correct==0, maxhdg=3;
    else, maxhdg = length(uhdg);
    end
    
    clear g L
    % (loop over h to allow animating figures one heading at a time)
    for h = 1:maxhdg      
%         g(h) = plot(squeeze(X(c,h,:)),squeeze(Y(c,h,:)),'color',cmap(h,:),'LineWidth', 2); hold on;
        g(h) = errorbar(squeeze(X(c,h,:)),squeeze(Y(c,h,:)),squeeze(Ye(c,h,:)),'color',cmap(h,:),'LineWidth', 2); hold on;

        set(g(h),'MarkerSize',10,'MarkerFaceColor',cmap(h,:));
        xlim([0.3 1.5]);
        
        if correct == -1, ylim([.35 .95]); % all trials
        else, ylim([0 1]); 
        end
           
        % set(gca,'Xtick',-2:1:2,'Ytick',-2:1:2);
        if c<=size(ucond,1) && ucond(c,2)==3,xlabel('RT (s)'); end 
        if mod(c,2)==0
            if conftask==1, ylabel('SEP'); else, ylabel('P(High Bet)'); end
        end
        changeAxesFontSize(gca,16,16); set(gca,'box','off');
        title(titles{c});
    end
    sh=suptitle('Confidence-RT'); set(sh,'fontsize',16,'fontweight','bold');
    
    if correct<0
        figure(20); 
        set(gcf,'Color',[1 1 1],'Position',[600 200 270*2 170*3],'PaperPositionMode','auto');
        subplot(3,2,subplotInd(c));
        
        clear g L
        % (loop over h to allow animating figures one heading at a time)
        for h = 1:length(uhdg)
            g(h) = errorbar(squeeze(X(c,h,:)),squeeze(Yc(c,h,:)),squeeze(Yce(c,h,:)),'color',cmap(h,:),'LineWidth', 2); hold on;
            
            set(g(h),'MarkerSize',10,'MarkerFaceColor',cmap(h,:));
            xlim([0.3 1.5]);
            ylim([0.35 1])
            
            % set(gca,'Xtick',-2:1:2,'Ytick',-2:1:2);
            if c<=size(ucond,1) && ucond(c,2)==3,xlabel('RT (s)'); end 
            if mod(c,2)==0, ylabel('P(correct)'); end
            changeAxesFontSize(gca,16,16); set(gca,'box','off');
            title(titles{c});
        end
        sh=suptitle('Accuracy-RT'); set(sh,'fontsize',16,'fontweight','bold');
    end
end




%{
nbins = 4; % how many RT quantiles?
inds  = floor(linspace(1,length(data.RT),nbins+1));
sortedRT = sort(data.RT);

edges = sortedRT(inds);
bins = discretize(data.RT,edges);

edges = [0 quantile(data.RT,nbins-1) Inf];
for q = 1:length(rtQ)-1
    J = theseRT>=rtQ(q) & theseRT<rtQ(q+1);
    xRT = mean(theseRT(J));
    Y(c,h,q) = mean(theseConf(J));
    Ye(c,h,q) = std(theseConf(J)) / sqrt(sum(J));
end

xRT = edges(1:end-1)+diff(edges)/2;

ahdgs = unique(abs(hdgs));
n = nan(nbins,length(mods),length(cohs),length(deltas),length(ahdgs));
Pcorrect = n;
PhighbetALL = n;
PhighbetCOR = n;
PhighbetERR = n;


for m=1:length(mods)
for c=1:length(cohs)
for d=1:length(deltas)
    for h=1:length(ahdgs)
        
        I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(d) & abs(data.heading)==ahdgs(h);
        
        for r=1:nbins
            R = bins==r;
            
            n(r,m,c,d,h) = sum(I & R);
            nC(r,m,c,d,h) = sum(I & R & data.correct);
            nE(r,m,c,d,h) = sum(I & R & ~data.correct);
            
            Pcorrect(r,m,c,d,h) = sum(I & R & data.correct) ./ n(r,m,c,d,h);
            
            if conftask==1
                PhighbetALL(r,m,c,d,h) = mean(data.conf(I & R));
                PhighbetCOR(r,m,c,d,h) = mean(data.conf(I & R & data.correct));
                PhighbetERR(r,m,c,d,h) = mean(data.conf(I & R & ~data.correct));
                
                

            elseif conftask==2 % PDW
            
                PhighbetALL(r,m,c,d,h) = sum(I & R & data.PDW==1) ./ n(r,m,c,d,h);
                PhighbetCOR(r,m,c,d,h) = sum(I & R & data.PDW==1 & data.correct) ./ nC(r,m,c,d,h);
                PhighbetERR(r,m,c,d,h) = sum(I & R & data.PDW==1 & ~data.correct) ./ nE(r,m,c,d,h);
            end
        end
    end
end
end
end

modlabels = {'Ves','Vis','Comb'};
mcols = {'Greys','Reds','Blues'};
D = 2;
figure('color','w')
for c=1:length(cohs)
    for m=1:length(mods)
        
        cmap = flipud(cbrewer('seq',mcols{m},length(hdgs)*2));
        cmap = cmap([9 7 5 3 1],:);

        subplot(length(mods),length(cohs),c+(m-1)*length(cohs)); hold on
        if m==1 && c~=1, delete(gca); continue, end

        temp = squeeze(Pcorrect(:,m,c,D,:));
        for h=1:size(temp,2)
        	plot(xRT,temp(:,h),'color',cmap(h,:),'marker','o','linew',2);  
        end
        
        if m~=1
            title([modlabels{m} ', coh = ' num2str(cohs(c))]); 
            if m==2 && c==1
                ylabel('P(Correct)');
            end
            if m==3
                xlabel('RT (s)')
            end
        else
            title(modlabels{m})
        end
        axis([min(xRT)-0.1 max(xRT)+0.1 0 1])
        changeAxesFontSize(gca,15,15);

    end
end


%%
temp1 = PhighbetCOR;
figure('color','w')
for c=1:length(cohs)
    for m=1:length(mods)
        
        cmap = flipud(cbrewer('seq',mcols{m},length(hdgs)*2));
        cmap = cmap([9 7 5 3 1],:);

        subplot(length(mods),length(cohs),c+(m-1)*length(cohs)); hold on
        if m==1 && c~=1, delete(gca); continue, end
        temp = squeeze(temp1(:,m,c,D,:));
        for h=1:size(temp,2)
        	plot(xRT,temp(:,h),'color',cmap(h,:),'marker','o','linew',2);
        end
        if m~=1
            title([modlabels{m} ', coh = ' num2str(cohs(c))]); 
            if m==2 && c==1
                if conftask==2
                    ylabel('P(High Bet)');
                elseif conftask==1
                    ylabel('Sacc EP');
                end
            end
            if m==3
                xlabel('RT (s)')
            end
        else
            title(modlabels{m})
        end
        axis([min(xRT)-0.1 max(xRT)+0.1 -0.3 1.3])
        changeAxesFontSize(gca,15,15);
    end
end
%}