%% getDataFromServer
% requires password(s) to be entered if a auth key not available 
% see e.g.:
% https://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/

% VERSION 2.0: 08-30-19 CF
% must be run by PLDAPS_preprocessing

% 10-2020: CF updated for OEphys and newline
% 08-2021: SJ updated for nexonar

isoephys = strcmp(remoteDir(end-5:end-1),'neuro');
isnexonar = strcmp(remoteDir(end-7:end-1),'nexonar'); 

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
if ~useVPN %any(strfind(ifstuff,'172.'))
    cmd = ['ssh fetschlab@172.30.3.33 ls ' remoteDir]; % MBI machine
else
    cmd = ['ssh fetschlab@10.161.240.133 ls ' remoteDir]; % probably off campus, try proxy IP (requires VPN)
end
% [~,remoteFileList] = system(cmd, '-echo');
[~,remoteFileList] = system(cmd);
if any(strfind(remoteFileList,'timed out')); error(remoteFileList); end

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
    
    if any(strfind(dateStr,thisDate)) && contains(remoteFiles{n},paradigm)
        currentFileList{m,1} = remoteFiles{n}; m=m+1;
        if ~any(strfind(localFileList,remoteFiles{n})) || overwriteLocalFiles % always copy if overwrite option selected
            if ~useVPN %any(strfind(ifstuff,'172.'))
                cmd = ['scp -r fetschlab@172.30.3.33:' remoteDir remoteFiles{n} ' ' localDir]; % MBI machine
            else
                cmd = ['scp -r fetschlab@10.161.240.133:' remoteDir remoteFiles{n} ' ' localDir]; % probably off campus, try proxy IP (requires VPN)
            end
            system(cmd,'-echo');
            
            % now clean it up and re-save
            if ~isoephys && ~isnexonar
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