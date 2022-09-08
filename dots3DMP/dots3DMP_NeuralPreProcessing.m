% Fetsch Lab
% S. J. Jerjian
% Updated August 2022
%
% this script creates a big datastruct containing spike sorted units and task data from rigB (dots3DMP recordings)
% user should specify the subject, dateRange, and paradigms
%
% organization of dataStruct is as follows:
%
% 1 - each row in dataStruct pertains to one recording 'session' (block/set) (in practice, this means each day)
% 2 - within a row, 
%       date and set
%       info field, with header information, 
%       data field
% data field is further split into fields for each paradigm recorded in that session. Each of these fields contains:
%
%   'events'    condition information for each trial, as well as the timing of key stimulus and behavioral events)
%   'pldaps'    useful pldaps information for each trial, largely bookkeeping
%   'units'     
%       'cluster_id'   :  1 x N ids, where N = number of clusters
%       'cluster_type' :  1 x N corresponding to cluster_ids above. By Kilosort convention,
%                               1 = MU, 2 = SU, (3 = noise, not included by default)
%                               if keepMU==1, MUs will be included, otherwise only SUs will be included in output struct
%       'spiketimes'   :  1 x N cell array containing spike times for each cluster, within the range of the given paradigm
%       'moreInfo'     :  1 x N cell array containing additional useful metadata      
% 
% **** A few notes for the user ****
% 
% Each unique row of dataStruct pertains to an individual 'set' of
% recordings, i.e. data recorded at a particular location. 
% % A 'good' recording day will likely only have one set, but 
% % in principle, there could be more than one set for a given day.
% Recordings are grouped by set, so multiple dots3DMP paradigms recorded in
% one set will be combined, whether the actual recording was done in
% multiple Trellis files or a single one - the "rec_group" group field in info struct is therefore critical!
% The timestamps for concatenated recordings are shifted according to the
% length of the overall data so that the range of events and spikes is
% matched for a given recording (nsEvents.analogData.timeStamps).
% i.e. if a recording is the first in the set, it's timestamps should start from around 0, whereas 
% % if it is later in the set they will start from some other time >0. The
% spiketimes will be matched accordingly, so the range of spiketimes should
% approximately match the range of e.g. nsEvents.Events.trStart.
%
% SJ 08-2022 added in metadata (getUnitInfo.m)
% SJ 08-2022 cleanUp option - to sub-select desirable units
% SJ 06-2022 significant updates
%            Fixed issues with processing of mksort data, 
%            + shifting of timestamps of multiple recordings. 
%            Switched dataCell from {} to dataStruct () struct format.
%
% Future improvements: 
% - include metadata about recorded units - depth, ch etc. (WORK IN
% PROGRESS 1) need to standardize for non-kilosort recordings too 2) need to store MDI and stopper depth)
% 
% - add option to load and append to existing dataStruct instead of
% generating from scratch at each run
% - consider whether any other info merits its own field in top level of
% dataStruct - MDI and stopper depth? distance to first contact, 'actual'
% distances between contacts
% could put this into the kilosort chanmap directly...instead of the
% unitless 1:32 currently used

% clear all
close all
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))

%% decide which files to load

% which paradigms do we care about?
paradigms = {'dots3DMPtuning','dots3DMP','RFMapping','VesMapping'};

% at the end, will discard any units whose gross spike rate is <minRate spikes
% per second, or the number of trials of any of the unique stimulus

runCleanUp = 0;
% conditions in that set is less than 5
parSelect  = {'dots3DMPtuning','dots3DMP'}; 
minRate    = 5;
minTrs     = 5;

subject = 'lucio';

% dateRange = [20220223:20220331 20220512:20220531];
% dateRange = 20220505:20220819;
dateRange = 20220901:20220902;

dateStr = num2str(dateRange(1));
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

%%
keepMU = 1;           % include all SU and MU
useSCP = 1;
useVPN = 1;
overwriteLocalFiles = 0; % set to 1 to always use the server copy

% SJ 04-2022
% download associated PDS data files (we'll need this for some cleanup)
% actually, just specify the mountDir, don't bother downloading the files

% localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/']; 
% remoteDir = ['/var/services/homes/fetschlab/data/' subject '/'];
% getDataFromServer;

mountDir  = ['/Volumes/homes/fetschlab/data/' subject '/'];    
PDSdir = mountDir; % reassign for later

localDir = ['/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/' subject '_neuro/'];
remoteDir = ['/var/services/homes/fetschlab/data/' subject '/' subject '_neuro/'];
%     mountDir  = ['/Volumes/homes/fetschlab/data/' subject '/' subject '_neuro/'];

%% grab the task events and info files
% this time we will locally download these (_RippleEvents and info files)

getNeuralEventsInfo; 

%% create the dataStruct

createSessionData;   

%% exclude cells which were not adequately recorded in ALL fundamental experiments
% parSelect

% let's remove units that don't have, say, at least 5 trials per unique
% condition?

if runCleanUp

fprintf('Cleaning up data, resaving\n')
dataStruct_clean = dataStruct;
removeEntireSession = false(size(dataStruct));

for s = 1:length(dataStruct)
    clear parUnits numSpikes numTrials enoughTrials

    % sessions that don't have all pars in parSelect get marked for
    % removal, but do it at the end so as not to mess up the loop counter
    if ~all(isfield(dataStruct(s).data,parSelect))
        removeEntireSession(s) = true;
        continue
    end

    for par = 1:length(parSelect)

        units  = dataStruct(s).data.(parSelect{par}).units;
        events = dataStruct(s).data.(parSelect{par}).events;

        parUnits(par,:)  = units.cluster_id;
        numSpikes(par,:) = cellfun(@length,units.spiketimes);
        numTrials(par,:) = length(events.trStart);


        stimCondList = [events.heading; events.modality; events.coherence; events.delta]';

        for u = 1:length(units.cluster_id)

            if ~isempty(units.spiketimes{u})
                [~,t] = min(abs(events.trStart-units.spiketimes{u}(1)));
                if ~isempty(t), itr_start=t; end
                [~,t] = min(abs(events.trStart-units.spiketimes{u}(end)));
                if ~isempty(t), itr_end=t; end
            end

            [uStimConds,~,ic]    = unique(stimCondList(itr_start:itr_end,:),'rows');
            [nTrConds,~]         = hist(ic,unique(ic));
            enoughTrials(par,u)  = all(nTrConds>=minTrs);
        end

    end

    numTrials       = repmat(numTrials,1,size(parUnits,2));

    parSpikeRate = numSpikes ./ numTrials;

    removeThese = any(parSpikeRate<minRate | ~enoughTrials,1);

    for par = 1:length(parSelect)
        
        units  = dataStruct(s).data.(parSelect{par}).units;
        units.cluster_id(removeThese) = [];
        units.cluster_type(removeThese) = [];
        units.spiketimes(removeThese) = [];

        % overwrite
        dataStruct_clean(s).data.(parSelect{par}).units = units;
    end
end

dataStruct_clean(removeEntireSession) = [];
dataStruct = dataStruct_clean;

file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '_neuralData_clean.mat'];
save([localDir(1:length(localDir)-length(subject)-7) file], 'dataStruct');

end