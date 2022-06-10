% last updated June 2022
% SJ

% this creates a big cell containing neural activity and information from
% rigB
% specify the subject, dateRange, and paradigms to draw from

% organization of dataCell is as follows:

% 1 - each cell in dataCell pertains to one recording 'session' (block/set) (in practice, this means each day)
% 2 - within a cell, 
%       info field, with header information, 
%       data field
% data field is further split into fields for each paradigm, with each
% paradigm field containing:
%   'events' condition information for each trial, as well as the timing of key stimulus and behavioral events)
%   'pldaps' useful pldaps information for each trial, largely bookkeeping
%   'cluster_id' 1 x N ids from Kilosort, where N = number of clusters
%   'cluster_type' 1 x N, 1 - MU, 2 - SU, 3 - noise (3 is skipped out, MU
%   will also be skipped if keepMU below == 0)
%   'spiketimes' 1 x N cell array containing spike times in seconds for
%   each cluster, within the range of the given paradigm

% example
% select the motionStartTimes for dots3DMP task in session 1, and
% spiketimes for first neuron

% s=1; u=1;
% motionStart = dataCell{s}.data.dots3DMP.events.stimOn;
% spktimes    = dataCell{s}.data.dots3DMP.spiketimes{u};



% TODO 
% - add option to append to existing dataCell
% - check that info contains all useful information
% - exclusion criteria for cells
%   - cluster_type, but also numSpikes/spikerate, presence in all paradigms
%
%

clear all
close all
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))

%% decide which files to load

% which paradigms do we care about
paradigms = {'dots3DMPtuning','dots3DMP','RFMapping','VesMapping'};
% paradigms = {'RFMapping'};

subject = 'lucio';

% dateRange = 20220218;
dateRange = 20220512:20220531;

dateStr = num2str(dateRange(1));
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

%%
useSCP = 1;
useVPN = 0;
overwriteLocalFiles = 1; % set to 1 to always use the server copy

% SJ 04-2022
% download associated PDS data files (we'll need this for some cleanup)
% actually, just specify the mountDir, don't bother downloading the files
% as my local mac is filling up with these files way too quickly
% localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/']; 
% remoteDir = ['/var/services/homes/fetschlab/data/' subject '/'];
mountDir  = ['/Volumes/homes/fetschlab/data/' subject '/'];
% getDataFromServer;
    
PDSdir = mountDir; % reassign for later

localDir = ['/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/' subject '_neuro/'];
remoteDir = ['/var/services/homes/fetschlab/data/' subject '/' subject '_neuro/'];
%     mountDir  = ['/Volumes/homes/fetschlab/data/' subject '/' subject '_neuro/'];

getNeuralEventsInfo; % grab the task events in rec time, and info files
createSessionData;   % create the dataCell
