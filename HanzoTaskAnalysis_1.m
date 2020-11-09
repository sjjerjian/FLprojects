% Preliminary analysis of Hanzo neural recordings
% SJ started 09/2020

% ch30 (31 on OpenEphys because chnums start at 0) on 09-24 had a unit sorted online
% recorded during RFmapping and dots task

sr = 30000;

%% load necessary data and clean up a bit

recfolder = '2020-10-06_15-19-44_exp';
% recfolder = '2020-09-24_14-29-21_e';

cd(sprintf('/Volumes/homes/fetschlab/data/hanzo_neuro/%s',recfolder));

% for ch=8:31
%     spikes_file = sprintf('SEp157.0n%d.spikes',ch);
%     [~,~,info]=load_open_ephys_data(spikes_file);
%     if any(info.sortedId)
%         ch
%     end
% end
    
ch=9;
spikes_file = sprintf('SEp157.0n%d.spikes',ch);    
[spikes,tSpike,info]=load_open_ephys_data(spikes_file);
[~, ~, ~, events] = load_open_ephys_data('messages.events', 1); 
[ttlBits, ttlTime] = load_open_ephys_data('all_channels.events');

%%
% load PLDAPS
% dont use clean version, as it is excluding trials missing choices which is causing
% mismatch in number of trials
% file = 'hanzo_20200924-20200925.mat';
% load(sprintf('Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/%s',file))
% trs = strcmp(data.filename,'hanzo20200924Dots1430');
% fnames = fieldnames(data);
% for F = 1:length(fnames)
%     data.(fnames{F}) = data.(fnames{F})(trs);
% end

file = 'hanzo20200924Dots1430.PDS';
load(sprintf('/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/hanzo/%s',file),'-mat') % load 'partly' cleaned-up version from local

% events vars are 1290 long, PDS is 1288 trials (why?)
% need to remove first two in events to match with PDS

fnames = fieldnames(events);
for F = 1:length(fnames)
    if strcmp(fnames{F},'TDirection') % why is this one shorter by 1?
        events.(fnames{F})(1) = [];
    else
        events.(fnames{F})(1:2) = [];
    end
end

%% 1. Plot waveforms of unit during task

% first 20 spikes
% last 20 spikes
% mean

spikes_sel = spikes(info.sortedId==3,:)./info.gain(info.sortedId==3)*1e6;
t = (0:size(spikes_sel,2)-1)/sr*1000;
figure('color','white'); hold on
plot(t,spikes_sel(1:20,:)','k','linew',0.25)
plot(t,spikes_sel(end-20:end,:),'color',[1 1 1]*0.5,'linew',0.25)
plot(t,mean(spikes_sel,1),'color','r','linew',2)

xlabel('time (ms)')
ylabel('V') % divided by gain above to get true amplitude, multiply by 1e6 to plot in uV, but not sure the amplitude is reasonable?

%%
nn = nan(length(PDS.conditions),1);
PDScoh = nn;
direction = nn;
choice = nn;
RT = nn;
PDW = nn;
correct = nn;
goodtrial = zeros(length(PDS.conditions),1);

fixTime = nn;
targTime = nn;

% extract useful variables to make code simpler
for i=1:length(PDS.conditions)
    PDScoh(i) = PDS.conditions{i}.stimulus.coherence;
    direction(i) = PDS.conditions{i}.stimulus.direction;
    
    if isfield(PDS.data{i}.behavior,'choice')
        choice(i) = PDS.data{i}.behavior.choice;
        RT(i) = PDS.data{i}.behavior.RT;
        PDW(i) = PDS.data{i}.behavior.PDW;
        correct(i) = PDS.data{i}.behavior.correct;
        goodtrial(i) = PDS.data{i}.behavior.goodtrial;
        
        % fixation time relative to dots onset
        fixTime(i) = PDS.data{i}.stimulus.timeFpEntered - PDS.data{i}.stimulus.timeDotsOn;
        
        % target onset relative to dots onset
        targTime(i) = -PDS.data{i}.stimulus.delayToDots;
    end
    
end

% signed coh
coh = PDScoh;
coh(coh<0.001) = 0;
coh(direction==180) = -coh(direction==180);

%% calculate average spike rate in motion interval (Dots on to RT)

% times of motion from ttlBits
ttlBitsIndex = find(ttlBits);
TimeDotsOn = ttlTime(ttlBitsIndex(1:2:end));
TimeDotsOff = ttlTime(ttlBitsIndex(2:2:end));

temp = events.MotionStart;
temp(~isnan(temp)) = TimeDotsOn;
TimeDotsOn = temp;

temp = events.MotionStart;
temp(~isnan(temp)) = TimeDotsOff;
TimeDotsOff = temp;
Duration = TimeDotsOff - TimeDotsOn; % should be the same as RT from PLDAPS...

%TimeDotsOn = events.MotionStart/sr;

TimeFix   = TimeDotsOn + fixTime;   % absolute time of fixation onset
TimeTargs = TimeDotsOn + targTime;  % absolute time of target onset
TimeRT    = TimeDotsOn + RT;        % absolute time of RT

spkcount = nan(size(TimeDotsOn));
spktimes = tSpike(info.sortedId==3);

for tr=1:length(spkcount)
    spkcount(tr) = sum(spktimes >= TimeDotsOn(tr) & spktimes <= TimeRT(tr));
end
avgspkrate = spkcount ./ RT;

% plot FR as function of signed coherence, independent of choice

ucohs = unique(coh);
cohFRs_mean = nan(length(ucohs),1);
cohFRs_se  = nan(length(ucohs),1);
ntrs = nan(length(ucohs),1);
for c=1:length(ucohs)
   ntrs(c) = sum(coh==ucohs(c));
   cohFRs_mean(c) = nanmean(avgspkrate(coh==ucohs(c)));
   cohFRs_se(c)  = nanstd(avgspkrate(coh==ucohs(c)))/sqrt(ntrs(c));
end

figure(2); errorbar(ucohs,cohFRs_mean,cohFRs_se)
title('Firing rate during motion (Dots On - RT)')
xlabel('Signed coh');
ylabel('Avg firing rate, sp/s')
ylim([20 50])

%% plot depending on choice/correct?
% what about PDW?

%{
ucohs = unique(coh);
cohFRs_mean = nan(length(ucohs),2);
cohFRs_se  = nan(length(ucohs),2);
ntrs = nan(length(ucohs),2);
for c=1:length(ucohs)
    for co=0:1
        ntrs(c,co+1) = sum(coh==ucohs(c) & correct==co);
        cohFRs_mean(c,co+1) = nanmean(avgspkrate(coh==ucohs(c) & correct==co));
        cohFRs_se(c,co+1)  = nanstd(avgspkrate(coh==ucohs(c) & correct==co))/sqrt(ntrs(c));
    end
end

figure; 
hold on
errorbar(ucohs,cohFRs_mean(:,2),cohFRs_se(:,2),'b')
errorbar(ucohs,cohFRs_mean(:,1),cohFRs_se(:,1),'r')

xlabel('Signed coh');
ylabel('Avg firing rate, sp/s')
% ylim([20 50])
%}

%% now do time-resolved firing rate, aligned to Dots on and RT separately
     
% this defines how much time before and after each event we want
% e.g. from 400ms before dots onset to 700ms after. and -700ms to +1000ms
% around RT
tmin = [-1.0 -0.7];
tmax = [0.7 1.0];
alignEvent = {TimeDotsOn, TimeRT};

bin = 0.05; % bin size, 50ms
smbins = 3; % smoothing boxcar, 3 bins (i.e. 150ms)

if smbins/2==floor(smbins/2), smbins=smbins+1; end
smbins2=floor(smbins/2);

Ntr = length(TimeRT);

clear FRs Bin_edges FRcohs_mean FRcohs_se
for i=1:length(tmin)

     %define x-axis,
     % go forward and backward from each event separately to make sure event
     % is in the middle of two bins
     if smbins==1
        edges=0:-bin:tmin(i);
        if edges(end)~=tmin(i), edges(end+1)=edges(end)-bin; end
        edges1=0:bin:tmax(i);
        if edges1(end)~=tmax(i), edges1(end+1)=edges1(end)+bin; end
       
     else
        edges=0:-bin:tmin(i);
        if edges(end)~=tmin(i), edges(end+1)=edges(end)-bin; end
        edges(end+1:end+smbins2)=edges(end)-1*bin:-bin:edges(end)-smbins2*bin;
        
        edges1=0:bin:tmax(i);
        if edges1(end)~=tmax(i), edges1(end+1)=edges1(end)+bin; end
        edges1(end+1:end+smbins2)=edges1(end)+1*bin:+bin:edges1(end)+smbins2*bin;  
     end
     
        edges=[edges(end:-1:1) edges1(2:end)];
        ne=length(edges);
        fr = nan(Ntr,ne-1);
    
    for itr=1:Ntr
        if isnan(alignEvent{i}(itr)), continue; end
        
        be=alignEvent{i}(itr)+edges(1);
        en=alignEvent{i}(itr)+edges(end);
        
        inds=find((spktimes>be-smbins*bin)&(spktimes<en+smbins*bin)); %spike times for the trial
        
        %     nspk=length(inds);
        if isempty(inds), continue; end
        index1=spktimes(inds)-alignEvent{i}(itr);
        
        fr(itr,:) = histcounts(index1,edges);
        
        if smbins>1
            fr(itr,:) = smooth(fr(itr,:),smbins);
        end
    end
    
%     fr(isnan(fr)) = 0;
    fr(:,end+1) = NaN; % if concatenating cells, keep different alignments separated. also matches length of bin edges
    FRs{i} = fr/bin; % sp/s
    
    if smbins==1
        Bin_edges{i} = edges;
    elseif smbins>1
        % shift bin x-values to centers if we've done smoothing
        Bin_edges{i} = edges+diff(edges(1:2))/2;
    end   
end

for c=1:length(ucohs)
    for i=1:length(FRs)
        Ntrs(c) = sum(coh==ucohs(c) & goodtrial);
        FRcohs_mean{i}(:,c) = nanmean(FRs{i}(coh==ucohs(c),:));
        FRcohs_se{i}(:,c)   = nanstd(FRs{i}(coh==ucohs(c),:))/sqrt(Ntrs(c));
    end
end

meanEvTimes = [nanmean(fixTime) 0 nanmean(RT)];

figure('position',[200 200 750 250],'color','w'); hold on
cmap = flipud(cbrewer('div','RdYlBu',length(ucohs)));

for i=1:length(FRs)
    
    % plot on separate axes

    subplot(1,2,i); hold on
    for c=1:length(ucohs)
        plot(Bin_edges{i},FRcohs_mean{i}(:,c),'color',cmap(c,:),'linew',1.25)
    end
    if i==1
        alignEvStr = 'Dots Onset';
        plot([0 0],[0 100],'k-'); text(0, 72, 'Dots','horizo','center');
        plot([nanmean(RT) nanmean(RT)],[0 100],'k--'); text(nanmean(RT), 72, 'RT','horizo','center');
        plot([nanmean(fixTime) nanmean(fixTime)],[0 100],'k--'); text(nanmean(fixTime), 72, 'Fixation','horizo','center');
        plot([nanmean(targTime) nanmean(targTime)],[0 100],'k--'); text(nanmean(targTime), 72, 'Targets','horizo','center');

    else 
        alignEvStr = 'Saccade Onset'; 
        
        plot([0 0],[0 100],'k-'); text(0, 72, 'RT','horizo','center');
        plot(-[nanmean(RT) nanmean(RT)],[0 100],'k--'); text(-nanmean(RT), 72, 'Dots','horizo','center');
        
        legend(cellstr(num2str(ucohs)))   
    end
    xlabel(sprintf('Time aligned to %s (s)',alignEvStr))
    ylabel('Firing rate (spikes/s)')
    axis([tmin(i)-0.05 tmax(i)+0.05 0 70])
    
    % or same axes (make tmin/tmax shorter to avoid excessive overlap)
    %{
    for c=1:length(ucohs)
        plot(Bin_edges{i}+meanEvTimes(i),FRcohs_mean{i}(:,c),'color',cmap(c,:),'linew',1.25)
    end
    plot([0 0],[0 100],'k-'); text(0, 62, 'Dots On','horizo','center');
    plot([nanmean(RT) nanmean(RT)],[0 100],'k-'); text(nanmean(RT), 62, 'RT','horizo','center');
     xlabel(sprintf('Time aligned to %s (s)',alignEvStr))
    ylabel('Firing rate (spikes/s)')
    ylim([0 60])
    %}
end

%% for given coherence or lumped cohs, plot FR as function of choice or PDW

% 1. absolute firing rate
% 2. residuals after mean subtraction 

% select coh
cc = 0; 
trs = abs(coh)==cc; % one coherence
% trs = abs(coh)<0.3 & abs(coh)>0.01;
%trs = true(size(coh));

cond = 'PDW'; % 'choice','PDW'
mean_subtract = 1;

clear FRcoh_cond
for i=1:length(FRs)
    
    switch cond
        case 'choice'
            
            FRcoh_cond{i}(:,1) = nanmean(FRs{i}(trs & choice==1,:));
            FRcoh_cond{i}(:,2) = nanmean(FRs{i}(trs & choice==2,:));
              
            if mean_subtract
                FRcoh_cond{i}(:,1) = nanmean(FRs{i}(trs & choice==1,:) - nanmean(FRs{i}(trs,:)));
                FRcoh_cond{i}(:,2) = nanmean(FRs{i}(trs & choice==2,:) - nanmean(FRs{i}(trs,:)));
            end
            
        case 'PDW'
            
            FRcoh_cond{i}(:,1) = nanmean(FRs{i}(trs & PDW==0,:));
            FRcoh_cond{i}(:,2) = nanmean(FRs{i}(trs & PDW==1,:));
            
            if mean_subtract
                FRcoh_cond{i}(:,1) = nanmean(FRs{i}(trs & PDW==0,:) - nanmean(FRs{i}(trs,:)));
                FRcoh_cond{i}(:,2) = nanmean(FRs{i}(trs & PDW==1,:) - nanmean(FRs{i}(trs,:)));
            end
    end
end

if mean_subtract
    yl = [-20 20];
else
    yl = [0 70];
end


figure('position',[200 200 750 250],'color','w'); 
cols = 'br';
for i=1:length(FRs)
    
    subplot(1,2,i); hold on
    for c=1:size(FRcoh_cond{i},2)
        plot(Bin_edges{i},FRcoh_cond{i}(:,c),'color',cols(c),'linew',1.25)
    end
    if i==1
        alignEvStr = 'Dots Onset';
        plot([0 0],[-100 100],'k-'); text(0, yl(2)+2, 'Dots','horizo','center');
        plot([nanmean(RT) nanmean(RT)],[-100 100],'k--'); text(nanmean(RT), yl(2)+2, 'RT','horizo','center');
        plot([nanmean(fixTime) nanmean(fixTime)],[-100 100],'k--'); text(nanmean(fixTime), yl(2)+2, 'Fixation','horizo','center');
        plot([nanmean(targTime) nanmean(targTime)],[-100 100],'k--'); text(nanmean(targTime), yl(2)+2, 'Targets','horizo','center');

    else 
        alignEvStr = 'Saccade Onset'; 
        
        plot([0 0],[-100 100],'k-'); text(0, yl(2)+2, 'RT','horizo','center');
        plot(-[nanmean(RT) nanmean(RT)],[-100 100],'k--'); text(-nanmean(RT), yl(2)+2, 'Dots On','horizo','center');
        
        switch cond
            case 'choice'
                text(0.2,yl(2)-5,sprintf('Leftward Choice, n = %d',sum(trs & choice==1)),'color',cols(1));
                text(0.2,yl(2)-10,sprintf('Rightward Choice, n = %d',sum(trs & choice==2)),'color',cols(2));
                
            case 'PDW'
                text(0.2,yl(2)-5,sprintf('Low Bet, n = %d',sum(trs & PDW==0)),'color',cols(1));
                text(0.2,yl(2)-10,sprintf('High Bet, n = %d',sum(trs & PDW==1)),'color',cols(2));
        end        
    end
    xlabel(sprintf('Time aligned to %s (s)',alignEvStr))
    ylabel('Firing rate (spikes/s)')
    
    if mean_subtract
        axis([tmin(i)-0.05 tmax(i)+0.05 yl])
    else
        axis([tmin(i)-0.05 tmax(i)+0.05 yl])
    end
end

