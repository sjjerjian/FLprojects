%% getDataFromServer
% requires password(s) to be entered if a auth key not available 
% see e.g.:
% https://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/

% VERSION 2.0: 08-30-19 CF
% must be run by PLDAPS_preprocessing

% 10-2020: CF updated for OEphys and newline

isoephys = strcmp(remoteDir(end-5:end-1),'neuro');

% get file list from local dir
if ~exist(localDir,'file')
    mkdir(localDir);
end
allFiles = dir(localDir);
if length(allFiles)>2
    localFileList = allFiles(3).name; % skip 1+2, they are are "." and ".."
    for n = 4:length(allFiles)
        localFileList = [localFileList newline allFiles(n).name];
    end
else
    localFileList = [];
end

% get file list from remote dir
cmd = 'ifconfig en7 inet'; % check if on local MBI network or need VPN workaround
[~,ifstuff] = system(cmd);
if any(strfind(ifstuff,'172.'))
    cmd = ['ssh fetschlab@172.30.3.33 ls ' remoteDir]; % MBI machine
else
    cmd = ['ssh fetschlab@10.161.240.133 ls ' remoteDir]; % probably off campus, try proxy IP (requires VPN)
end
[~,remoteFileList] = system(cmd);
if any(strfind(remoteFileList,'timed out')); error(remoteFileList); end


% % % % this should be obsolete if using newline, but hang onto it for a bit just in case
% % 
% % % filenameStart = strfind(lower(remoteFileList),lower(subject)); % lower makes it case-insensitive
% % % filenameEnd = strfind(remoteFileList,'.')+3; % . instead of .PDS because sometimes they are .mat
% % % if length(filenameStart) ~= length(filenameEnd)
% % %     error('invalid file list: starts and ends don''t match');
% % %     % now this should only happen if filename has >1 or <1 dot, or there
% % %     % are files with the wrong subject name (not just a case mismatch)
% % % end
% % % 
% % % % check each file for a match to the desired date range and paradigm,
% % % % and if it doesn't already exist locally, copy it over
% % % remoteFiles = cell(length(filenameStart),1);
% % % for n = 1:length(filenameStart)
% % %     remoteFiles{n} = remoteFileList(filenameStart(n):filenameEnd(n));
% % %     dateStart = strfind(remoteFiles{n},'20');
% % %     thisDate = remoteFiles{n}(dateStart(1):dateStart(1)+7); % (1) just in case '20' appears in a timestamp
% % %     dot = strfind(remoteFiles{n},'.');
% % %     thisProt = remoteFiles{n}(dateStart(1)+8:dot-5);
% % %     if strcmp(paradigm,'Dots') % kluge for Hanzo, don't require an exclusive match
% % %         foundProt = contains(remoteFiles{n},thisProt);
% % %     else
% % %         foundProt = strcmp(paradigm,thisProt);
% % %     end
% % %     if any(strfind(dateStr,thisDate)) && foundProt
% % %         if ~any(strfind(localFileList,remoteFiles{n})) || overwriteLocalFiles % always copy if overwrite option selected
% % % %             cmd = ['scp fetschlab@172.30.3.33:' remoteDir remoteFiles{n} ' ' localDir];
% % %             cmd = ['scp fetschlab@10.161.240.133:' remoteDir remoteFiles{n} ' ' localDir]; % temp: VPN workaround
% % %             system(cmd,'-echo');
% % %             % now clean it up and re-save
% % %             pdsCleanup
% % %         else
% % %             disp([remoteFiles{n} ' exists locally, not copied']);
% % %         end
% % %     end
% % % 
% % % end

newlines = strfind(remoteFileList,newline);

% check each file for a match to the desired date range and paradigm,
% and if it doesn't already exist locally, copy it over
remoteFiles = cell(length(newlines),1);
m = 1;
for n = 1:length(newlines)
    if n==1
        remoteFiles{n} = remoteFileList(1:newlines(1)-1);
    else
        remoteFiles{n} = remoteFileList(newlines(n-1)+1:newlines(n)-1);
    end
    dateStart = strfind(remoteFiles{n},'20');
    if isempty(dateStart); continue; end
    if isoephys
        thisDate = remoteFiles{n}(dateStart(1):dateStart(1)+9); % dateStart(1) just in case '20' appears in a timestamp
        thisDate([5 8]) = []; % OpenEphys dates are always yyyy-mm-dd, so remove the dashes
    else
        thisDate = remoteFiles{n}(dateStart(1):dateStart(1)+7); % dateStart(1) just in case '20' appears in a timestamp
    end
    
    if any(strfind(dateStr,thisDate))
        currentFileList{m,1} = remoteFiles{n}; m=m+1;
        if ~any(strfind(localFileList,remoteFiles{n})) || overwriteLocalFiles % always copy if overwrite option selected
            if any(strfind(ifstuff,'172.'))
                cmd = ['scp -r fetschlab@172.30.3.33:' remoteDir remoteFiles{n} ' ' localDir]; % MBI machine
            else
                cmd = ['scp -r fetschlab@10.161.240.133:' remoteDir remoteFiles{n} ' ' localDir]; % probably off campus, try VPN workaround
            end
            system(cmd,'-echo');
            
            % now clean it up and re-save
            if ~isoephys
                pdsCleanup
            end
        else
            disp([remoteFiles{n} ' exists locally, not copied']);
            if isoephys
                warning('Partially copied subfolders on localDir will be skipped!');
            end
        end
        % clean up oephys even if not copied
        if isoephys
            oephysCleanup
        end

    end
end

disp('done.');