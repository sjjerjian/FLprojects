function nex = processNexonarUDP(data)

% clear;clc
% 
% hdr.dataPath = 'C:\Users\fetschlab\data\';
% hdr.nexFilename =  uigetfile([hdr.dataPath '*.txt']);
% data = load([hdr.dataPath hdr.nexFilename]); 

%% Extract trial information and trial data
trialSeparators = [0; find(all(data==-1,2))];

% unique_nums = data(trialSeparators(1:end-1)+1,2);
% [~,ia] = unique(unique_nums,'stable');
% trialSeparators = trialSeparators(ia);

ntrials = length(trialSeparators)-1;
nex = struct('nexdata',[],'pldaps',[],'conditions',[],'behavior',[]);

for n = 1:ntrials
    trial_data = data( trialSeparators(n) + 2 : trialSeparators(n+1)-2, 1:5);
    
%     badrows = find(all(trial_data==-1,2));
%     if ~isempty(badrows)
%         badrows = [-1 0 1] + badrows;
%         trial_data(badrows,:) = [];
%     end
%     nex.nexdata{n} = trial_data;
%     
    trial_info = data( trialSeparators(n) + 1, :);
    
    if n == ntrials
        trial_outcome = data(end-1,:);
    else
        trial_outcome = data( trialSeparators(n+1)-1, :);
    end
    
    if ~ismember(trial_outcome(2), [0 1])
        % sometimes the trial info for the first trial is not sent (why?)
        % and if PLDAPS experiment stops due to an error (rather than
        % intentional quit), the data stream for the last trial will be cut
        % off. so we may miss trial data for the occasional trial, almost
        % always the first or last one in a block
        fprintf(' trial %d/%d, outcome data missing or corrupted...\n',n,ntrials)
        nex.behavior.goodtrial(n)   = 0;
        nex.behavior.choice(n)      = NaN;
        nex.behavior.correct(n)     = NaN;
        nex.behavior.RT(n)          = NaN;
        nex.behavior.conf(n)        = NaN;
    else
        nex.behavior.goodtrial(n)   = trial_outcome(2);
        nex.behavior.choice(n)      = trial_outcome(3);
        nex.behavior.correct(n)     = trial_outcome(4);
        nex.behavior.RT(n)          = trial_outcome(5);
        nex.behavior.conf(n)        = trial_outcome(6);
    end
    
    nex.nexdata{n} = trial_data;
    nex.pldaps.iTrial(n) = trial_info(1);
    nex.pldaps.trialSeed(n) = trial_info(2);
    
    
    nex.conditions.modality(n)  = trial_info(3);
    nex.conditions.heading(n)   = trial_info(4);
    nex.conditions.coherence(n) = trial_info(5);
    nex.conditions.delta(n)     = trial_info(6);
    
    % for VesMapping, need to find a way to switch this automatically
    % (specify paradigm)
%     nex.conditions.headingFreq(n)  = trial_info(3);
%     nex.conditions.headingAmpl(n)   = trial_info(4);
%     nex.conditions.headingTheta(n) = trial_info(5);
%     nex.conditions.headingPhi(n)     = trial_info(6);
    
    
end

% some kluge cleanups
% for test20220922...why is 
% nex.nexdata{33} = [nex.nexdata{32}; nex.nexdata{33}];
% nex.nexdata(32) = [];

% fnames = fieldnames(nex.conditions);
% for f = 1:length(fnames)
%     nex.conditions.(fnames{f})(32) = [];
% end
% fnames = fieldnames(nex.pldaps);
% for f = 1:length(fnames)
%     nex.pldaps.(fnames{f})(32) = [];
% end
% fnames = fieldnames(nex.behavior);
% for f = 1:length(fnames)
%     nex.behavior.(fnames{f})(32) = [];
% end


% dot = strfind(hdr.nexFilename,'.');
% hdr.filename = [hdr.nexFilename(1:dot-1) '_nexonar.mat'];
% save([hdr.dataPath hdr.filename],'nex','hdr')
%% 

if 0
    
F = fieldnames(newdata);
for f = 1:length(F)
    newdata.(F{f})(goodtrials==0) = [];
end
nexdata(goodtrials==0) = [];

%% plots
XYZ = 'xyz';
cols = {'k','r','b'};
% ymax = [
for t=1:length(nexdata)
    for xyz=1:3
        subplot(2,3,xyz); hold on;
        plot(nexdata{t}(:,1),nexdata{t}(:,2+xyz)-nexdata{t}(1,2+xyz),'color',cols{newdata.modality(t)},'linew',1.5);
        xlim([0 6000])
        title(sprintf('%s displacement [cm]',XYZ(xyz)))
        if xyz==1
            ylim([-35 35])
        elseif xyz==2
            ylim([-10 110])
        else
            ylim([-10 130]);
        end
        
        subplot(2,3,xyz+3); hold on;
        plot(nexdata{t}(:,1),gradient(nexdata{t}(:,2+xyz)-nexdata{t}(1,2+xyz)),'color',cols{newdata.modality(t)},'linew',1.5);
        xlim([0 6000])
    end
end


%%
mods  = unique(newdata.modality);
cohs  = unique(newdata.coherence);
deltas = unique(newdata.delta);
hdgs  = unique(newdata.heading);

maxlen = max(cellfun(@length,nexdata));
for m=1:length(mods)
    for h=1:length(hdgs)
        temp = nexdata(newdata.heading==hdgs(h) & newdata.modality==mods(m));
        temptrs = nan(maxlen,length(temp), 3);
        for tr = 1:length(temp)
            temptrs(1:size(temp{tr},1), tr, 1:3) = temp{tr}(:,3:5);
        end
        meanNexTraj(:,:,m,h) = squeeze(nanmean(temptrs,2));
    end
end

meanNexTraj = meanNexTraj - meanNexTraj(1,:,:,:);

ss = 2/maxlen;
figure;
cols = 'krb';
for xyz=1:3
    subplot(3,1,xyz); hold on;
    for m=1:3
        plot(squeeze(meanNexTraj(:,xyz,m,:)),'color',cols(m))
    end
end
    
end


