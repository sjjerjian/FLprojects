

filename = 'test20220120dots3DMP1730';

load([filename '_nexonar.mat']);
load([filename '.mat']);

%%

data = createDataStructure_oneFile(PDS,'lucio');
[nexPDS,nexClean,exitflag] = dots3DMP_nexonarCleanUp(nex,PDS);

data.nexonar = nexPDS';

%%
%%%%%% PLOTTING
mods = [1 2 3];
% hdgs = [-12 -6 -3 -1.5 0 1.5 3 6 12];
hdgs = unique(data.heading);

hdgCol = flipud(cbrewer('div','RdBu',length(hdgs)));

% set small hdg to zero, change color of zero heading to mark it easily
data.heading(abs(data.heading)<0.01) = 0;
hdgs(abs(hdgs)<0.01) = 0;
hdgCol(hdgs==0,:) = [0.5 0.5 0.5];

maxlen = max(cellfun(@(x) size(x,1), nex.nexdata));

% nexDat = nex.nexdata;
nexDat = data.nexonar;
nexMat = nan(maxlen,size(nexDat{1},2),length(nexDat));

for t=1:length(nexDat)
    nexMat(1:size(nexDat{t},1),:,t) = nexDat{t};
end

nexMat = permute(nexMat,[1 3 2]);
mu = nexMat(1,:,:);
nexMat = bsxfun(@minus,nexMat,mu);
% nexMat(:,:,3) = -nexMat(:,:,3);

nT = 200; %

modlabels = {'Ves','Vis','Comb'};
figure('position',[100 100 800 800],'color','w');
for m=1:length(mods)
    for h=1:length(hdgs)
        trs = ~isnan(data.choice) & data.heading==hdgs(h) & data.modality==mods(m) & data.delta==0;
        subplot(3,3,m); hold on; axis([-40 40 0 180]); offsetAxes; title(sprintf('%s,x-y',modlabels{m}));
        plot(nexMat(1:nT,trs,3),nexMat(1:nT,trs,4),'color',hdgCol(h,:),'linewidth',0.5);
        subplot(3,3,m+3); hold on; axis([-40 40 0 5]); offsetAxes; title(sprintf('%s,x-z',modlabels{m}));
        plot(nexMat(1:nT,trs,3),nexMat(1:nT,trs,5),'color',hdgCol(h,:),'linewidth',0.5);
        subplot(3,3,m+6); hold on; axis([0 180 -2 5]); offsetAxes; title(sprintf('%s,y-z',modlabels{m}));
        plot(nexMat(1:nT,trs,4),nexMat(1:nT,trs,5),'color',hdgCol(h,:),'linewidth',0.5);
    end
end

figure('position',[100 100 800 800],'color','w');
for m=1:length(mods)
    for h=1:length(hdgs)
        trs = ~isnan(data.choice) & data.heading==hdgs(h) & data.modality==mods(m) & data.delta==0;
        subplot(3,3,m); hold on; axis([0 5000 -40 40]); offsetAxes; title(sprintf('%s,T-x',modlabels{m}));
        plot(nexMat(:,trs,1),nexMat(:,trs,3),'color',hdgCol(h,:),'linewidth',0.5);
        subplot(3,3,m+3); hold on; axis([0 5000 0 160]); offsetAxes; title(sprintf('%s,T-y',modlabels{m}));
        plot(nexMat(:,trs,1),nexMat(:,trs,4),'color',hdgCol(h,:),'linewidth',0.5);
        subplot(3,3,m+6); hold on; axis([0 5000 -2 5]); offsetAxes; title(sprintf('%s,T-z',modlabels{m}));
        plot(nexMat(:,trs,1),nexMat(:,trs,5),'color',hdgCol(h,:),'linewidth',0.5);
    end
end

