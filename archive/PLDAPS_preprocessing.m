% grabs PLDAPS data from NAS as specified, then loads desired variables
% from .PDS files into a simpler data structure

% CF, adapted from offline analysis wrapper born 10-3-18

% VERSION 2.0: 08-30-19

% requires password(s) to be entered if a auth key not available 
% see e.g.:
% https://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/

clear all
close all

%% decide which files to load

% subject = 'hanzo';
% paradigm = 'Dots';

% dateRange = 20210208:20210212; % last week

% % month by month
% dateRange = 20190301:20190331; % skip, too early?
% dateRange = 20190401:20190430; % early days; some bias, no slope eff -- but got to V shape pretty fast!
% dateRange = 20190501:20190531; % fine
% dateRange = 20190601:20190630; % fine
% dateRange = 20190701:20190731; % fine (**first RT** (not much))
% dateRange = 20190801:20190831; % good
        % no data in september 2019
% dateRange = 20191001:20191031; % good PDW curve (no RT), but opposite bias as below 
% dateRange = 20191101:20191130; % fantastic! (no RT)
        % no data in december 2019
% dateRange = 20200101:20200131; % some bias
% dateRange = 20200201:20200229; % large bias
% dateRange = 20200301:20200331; % small bias but no slope effect
        % no data in april or may 2020
% dateRange = 20200601:20200630; % small bias, good slope effect, but flat PDW
% dateRange = 20200701:20200731; % bias is back, otherwise improving
% dateRange = 20200801:20200831; % good. interesting RT separated by conf!
% dateRange = 20200901:20200930; % good. that RT effect is now gone.
% dateRange = 20201001:20201031; % good. RT effect back (modest). PDW pretty shallow
% dateRange = 20201101:20201130; % beautiful
% dateRange = 20201201:20201231; % not as nice but small N

% dateRange = 20190401:20201231; % all 'good' data

% dateRange = 20190401:20191130; % 'early'
% dateRange = 20200801:20201231; % 'late';

% % Warning: error loading hanzo20191011Dots1341.PDS
% % Warning: error loading hanzo20200910Dots1221.PDS
% % Warning: error loading hanzo20201102Dots1409.PDS


% dateRange = [20190501:20190831 20191101:20191130];

% dateRange = 20200801:20201130; % good!




% now do Genji!


% subject = 'human';
% paradigm = 'dots3DMP';
% dateRange = 20190612:20191231; % pre-RT

subject = 'human';
paradigm = 'dots3DMP';
dateRange = 20200213:20210526; % RT

% subject = 'human';
% paradigm = 'dots3DMP';
% dateRange = 20190612:20200308; % everything!

% subject = 'lucio';
% paradigm = 'dots3DMP';
% dateRange = 20200315:20210512; % RT

dateStr = num2str(dateRange(1));
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

% localDir = ['/Users/chris/Documents/MATLAB/PLDAPS_data/' subject '/'];
localDir = ['/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/' subject '/'];
% 
% remoteDir = ['/var/services/homes/fetschlab/data/' subject '_basic/'];
remoteDir = ['/var/services/homes/fetschlab/data/' subject 'dots3DMP_basic/'];


%% get PDS files from server -- DON'T FORGET VPN
% will skip files that already exist locally, unless overwrite set to 1

useVPN = 0;
overwriteLocalFiles = 0; % set to 1 to always use the server copy
getDataFromServer % now also includes pdsCleanup to reduce file size and complexity


%% create data structure

createDataStructure


%% optional: save data struct to a mat file so you don't have to repeat the time consuming step

if sum(diff(dateRange)>1)==0
    file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
elseif sum(diff(dateRange)>1)==1
    file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(diff(dateRange)>1)) '+' num2str(dateRange(find(diff(dateRange)>1)+1)) '-' num2str(dateRange(end)) '.mat'];
else
    warning('Multiple breaks in date range, could overwrite existing file with different intervening dates!');
    disp('Press a key to continue');
    pause
    file = [subject '_' num2str(dateRange(1)) '---' num2str(dateRange(end)) '.mat'];
end
    


try
data = rmfield(data,'dotPos'); % CAREFUL
catch
end

disp('saving...');
save([localDir(1:length(localDir)-length(subject)-1) file], 'data');

% otherwise for larger files will need: 
% save([localDir(1:length(localDir)-length(subject)-1) file], 'data','-v7.3');

disp('done.');

