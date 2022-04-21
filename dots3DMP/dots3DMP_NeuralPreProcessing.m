
clear all
close all
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))

%% decide which files to load

paradigm = 'dots3DMPtuning';
today    = str2double(datestr(now,'yyyymmdd'));


subject = 'lucio';
% dateRange = 20220101:today; % RT
% maybe need to modify this to easily skip certain dates if desired...
dateRange = 20220218:20220308;

dateStr = num2str(dateRange(1));
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

%%
useSCP = 1;
useVPN = 0;
overwriteLocalFiles = 0; % set to 1 to always use the server copy

% download associated PDS data files (we'll need this for some cleanup)
localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/']; 
remoteDir = ['/var/services/homes/fetschlab/data/' subject '/'];
mountDir  = ['/Volumes/homes/fetschlab/data/' subject '/'];
% getDataFromServer;
    
PDSdir = mountDir; % reassign for later
%%
localDir = ['/Users/stevenjerjian/Desktop/FetschLab/Analysis/data/' subject '_neuro/'];
remoteDir = ['/var/services/homes/fetschlab/data/' subject '/' subject '_neuro/'];
%     mountDir  = ['/Volumes/homes/fetschlab/data/' subject '/' subject '_neuro/'];

% just the RippleEvents and info files
createNeuralDataStruct;
