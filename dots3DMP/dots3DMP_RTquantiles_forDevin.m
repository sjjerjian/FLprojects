load('lucio_20210401-20210510_clean.mat')

%%
% hack to remove some long RTs which are skewing the data, for now
removethese = data.RT>1.2;
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

nbins = 6; % how many RT quantiles?
inds  = floor(linspace(1,length(data.RT),nbins+1));

%%%% INSERT CODE HERE TO SORT THE RT DATA
sortedRT = ...

edges = sortedRT(inds);
bins = discretize(data.RT,edges);

xRT = edges(1:end-1)+diff(edges)/2;

ahdgs = unique(abs(hdgs));
n = nan(nbins,length(mods),length(cohs),length(deltas),length(ahdgs));

Pcorrect = n;
PhighbetALL = n;
PhighbetCOR = n;
PhighbetERR = n;

%% CALCULATE PROPORTION CORRECT AND PROPORTION HIGH BET for each condition
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
            PhighbetALL(r,m,c,d,h) = sum(I & R & data.PDW==1) ./ n(r,m,c,d,h);
            PhighbetCOR(r,m,c,d,h) = sum(I & R & data.PDW==1 & data.correct) ./ nC(r,m,c,d,h);
            PhighbetERR(r,m,c,d,h) = sum(I & R & data.PDW==1 & ~data.correct) ./ nE(r,m,c,d,h);

        end
    end
end
end
end

%% PLOTTING

modlabels = {'Ves','Vis','Comb'};
cols = {'krb'};
figure('color','w')
for c=1:length(cohs)
    for m=1:length(mods)

        subplot(length(mods),length(cohs),c+(m-1)*length(cohs)); hold on
        if m==1 && c~=1, delete(gca); continue, end
        for h=1:size(temp,2)
            %%%% INSERT CODE HERE to plot Pcorrect for each condition vs RT.
            
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
        axis([min(xRT)-0.1 max(xRT)+0.1 0.4 1])
        changeAxesFontSize(gca,15,15);

    end
end

%%%% INSERT CODE HERE to do a similar figure but for PHighBets...
