

subject = 'lucio';
paradigm = 'dots3DMP';
dateRange = 20210714;

dateStr = num2str(dateRange(1));
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/nexonar/'];
remoteDir = ['/var/services/homes/fetschlab/data/' subject '/' subject '_nexonar/'];
% remoteDir = ['/var/services/homes/fetschlab/data/'];

useVPN = 0;
overwriteLocalFiles = 1; % set to 1 to always use the server copy
getDataFromServer % now also includes pdsCleanup to reduce file size and complexity

%%
localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/'];
% remoteDir = ['/var/services/homes/fetschlab/data/' subject '/' ];

localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/'];
remoteDir = ['/var/services/homes/fetschlab/data/'];

getDataFromServer % now also includes pdsCleanup to reduce file size and complexity


%%
clear all

subject = 'lucio';
localDirNex = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/nexonar/'];
% localDirPDS = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/'];

filename = 'lucio20210714dots3DMP1252';
load(fullfile(localDirNex,[filename '_nexonar.mat']))
% % load(fullfile(localDirPDS,[filename '.mat']))

% clean-up bug in nex trial storage
% badtrialSeed = find(nex.pldaps.trialSeed==7);
%%
% kluge cleanup some old data for Arielle 07-19-2021

lens= cellfun(@(x) size(x,1),nex.nexdata);

wrong_iTrials = find(nex.pldaps.trialSeed==7);
nex.pldaps.iTrial(wrong_iTrials) = nex.pldaps.iTrial(wrong_iTrials-1)+1;

split_trials = find(diff(nex.pldaps.iTrial)==0);

% trialSeedSplit(:,1) = nex.pldaps.trialSeed(split_trials);
% trialSeedSplit(:,2) = nex.pldaps.trialSeed(split_trials+1);

discardInds = split_trials+1;
%%%%
for i=1:length(split_trials)
    nex.nexdata{split_trials(i)} = [nex.nexdata{split_trials(i)}; nex.nexdata{split_trials(i)+1}];
    nex.nexdata(split_trials(i)+1) = [];
end

fnames = fieldnames(nex.conditions);
for f=1:length(fnames)
    nex.conditions.(fnames{f})(split_trials+1) = [];
end

fnames = fieldnames(nex.behavior);
for f=1:length(fnames)
    nex.behavior.(fnames{f})(split_trials+1) = [];
end


fnames = fieldnames(nex.pldaps);
for f=1:length(fnames)
    nex.pldaps.(fnames{f})(split_trials+1) = [];
end

%%
clearvars -except nex hdr


