addpath(genpath('/Users/stevenjerjian/Desktop/PhD/Codes/General/'))
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/Codes/offlineTools'))

cd('/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data');

%%
% load('lucio_20210315-20210512_clean.mat')
load('human_20200213-20210512_clean.mat')

%%
% hack to remove some long RTs which are skewing the data, for now
removethese = data.RT>2.3;
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

nbins = 8; % how many RT quantiles?
inds  = floor(linspace(1,length(data.RT),nbins+1));
sortedRT = sort(data.RT);

edges = sortedRT(inds);
bins = discretize(data.RT,edges);

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
        cmap = cmap([7 5 3 1],:);

        subplot(length(mods),length(cohs),c+(m-1)*length(cohs)); hold on
        if m==1
            if c>1, delete(gca); continue, end
        end
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
temp1 = PhighbetERR;
figure('color','w')
for c=1:length(cohs)
    for m=1:length(mods)
        
        cmap = flipud(cbrewer('seq',mcols{m},length(hdgs)*2));
        cmap = cmap([7 5 3 1],:);

        subplot(length(mods),length(cohs),c+(m-1)*length(cohs)); hold on
        if m==1 && c~=1, delete(gca); continue, end
        temp = squeeze(temp1(:,m,c,D,:));
        for h=1:size(temp,2)
        	plot(xRT,temp(:,h),'color',cmap(h,:),'marker','o','linew',2);
        end
        if m~=1
            title([modlabels{m} ', coh = ' num2str(cohs(c))]); 
            if m==2 && c==1
%                 ylabel('P(High Bet)');
                ylabel('Sacc EP');

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