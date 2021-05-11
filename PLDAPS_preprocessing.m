% grabs PLDAPS data from NAS as specified, then loads desired variables
% from .PDS files into a simpler data structure

% CF, adapted from offline analysis wrapper born 10-3-18

% VERSION 2.0: 08-30-19

% requires password(s) to be entered if a auth key not available 
% see e.g.:
% https://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/

clear all

%% decide which files to load

subject = 'lucio';
paradigm = 'dots3DMP';
dateRange = 20210315:20210510; % 
% % 
% % % Warning: error loading hanzo20191011Dots1341.PDS

% subject = 'human';
% paradigm = 'dots3DMP';
% dateRange = 20190612:20191231; % pre-RT

% subject = 'human';
% paradigm = 'dots3DMP';
% dateRange = 20200213:20200308; % RT

%%
dateStr = num2str(dateRange(1));
if length(dateStr)>1
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end
end

% localDir = ['/Users/chris/Documents/MATLAB/PLDAPS_data/' subject '/'];
localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/'];
remoteDir = ['/var/services/homes/fetschlab/data/' subject '_basic/'];

useVPN = 0; % set to 1 if connected to JHU VPN
%createLocalFiles = 1; % set to 0 to run pdsCleanup on server
% server pdsCleanup not yet implemented...trying to solve error issue with
% 'load' and 'save' of PDS file on remote server

overwriteLocalFiles = 0; % set to 1 to always use the server original copy

%% get PDS files from server
% will skip files that already exist locally, unless overwrite set to 1

getDataFromServer % includes pdsCleanup to reduce file size and complexity

%% create data structure

createDataStructure


%% optional: save data struct to a mat file so you don't have to repeat the time consuming step
file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
save([localDir(1:length(localDir)-length(subject)-1) file], 'data');
