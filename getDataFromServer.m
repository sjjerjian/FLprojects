%% getDataFromServer
% requires password(s) to be entered if a auth key not available 
% see e.g.:
% https://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/

% VERSION 2.0: 08-30-19 CF
% must be run by PLDAPS_preprocessing


% get file list from local dir
allFiles = dir(localDir);
localFileList = allFiles(3).name; % skip 1+2, they are are "." and ".."
for n = 4:length(allFiles)
    localFileList = [localFileList newline allFiles(n).name];
end

% get file list from remote dir
cmd = ['ssh fetschlab@172.30.3.33 ls ' remoteDir];
[~,remoteFileList] = system(cmd);
filenameStart = strfind(remoteFileList,subject);
filenameEnd = strfind(remoteFileList,'.PDS')+3;
if length(filenameStart) ~= length(filenameEnd)
    error('invalid file list: starts and ends don''t match');
        % this can happen, for example, if not all the files in the folder
        % match the subject name, case sensitive (e.g. Hanzo vs. hanzo), or
        % if some of them are not .PDS
end

% check each file for a match to the desired date range and paradigm,
% and if it doesn't already exist locally, copy it over
remoteFiles = cell(length(filenameStart),1);
for n = 1:length(filenameStart)
    remoteFiles{n} = remoteFileList(filenameStart(n):filenameEnd(n));
    dateStart = strfind(remoteFiles{n},'20');
    thisDate = remoteFiles{n}(dateStart(1):dateStart(1)+7); % (1) just in case '20' appears in a timestamp
    dot = strfind(remoteFiles{n},'.');
    thisProt = remoteFiles{n}(dateStart(1)+8:dot-5);
    if strcmp(paradigm,'Dots') % kluge for Hanzo, don't require an exclusive match
        foundProt = contains(remoteFiles{n},thisProt);
    else
        foundProt = strcmp(paradigm,thisProt);
    end
    if any(strfind(dateStr,thisDate)) && foundProt
        if ~any(strfind(localFileList,remoteFiles{n})) || overwriteLocalFiles % always copy if overwrite option selected
            cmd = ['scp fetschlab@172.30.3.33:' remoteDir remoteFiles{n} ' ' localDir];
            system(cmd,'-echo');
            % now clean it up and re-save
            pdsCleanup
        else
            disp([remoteFiles{n} ' exists locally, not copied']);
        end
    end
    
end

disp('done.');