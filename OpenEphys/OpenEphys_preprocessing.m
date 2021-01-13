% grabs OpenEphys data from NAS as specified, then loads desired variables
% from .PDS files into a simpler data structure

% 10-18-20

% requires password(s) to be entered if a auth key not available 
% see e.g.:
% https://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/

clear all

%% decide which files to load

subject = 'hanzo';
paradigm = 'Dots';
dateRange = 20201201:20201202;

dateStr = num2str(dateRange(1));
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

localDir = ['/Users/chris/Documents/MATLAB/OpenEphys_data/' subject '/'];
remoteDir = ['/var/services/homes/fetschlab/data/' subject '_neuro/'];


%% get PDS files from server
% will skip files that already exist locally, unless overwrite set to 1

overwriteLocalFiles = 0; % set to 1 to always use the server copy
getDataFromServer % now also includes pdsCleanup to reduce file size and complexity


%% data browser (work in progress, may be abandoned)

oephysDataBrowser





% TBD:
%
% createDataStructure
%
% 
% %% optional: save data struct to a mat file so you don't have to repeat the time consuming step
% file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
% 
% data = rmfield(data,'dotPos'); % CAREFUL
% 
% save([localDir(1:length(localDir)-length(subject)-1) file], 'data');
% 
% % otherwise for larger files will need: 
% % save([localDir(1:length(localDir)-length(subject)-1) file], 'data','-v7.3');
% 
% 
% disp('done.');
% 
