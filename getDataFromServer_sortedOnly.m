%% getDataFromServer_sortedOnly
% requires password(s) to be entered if a auth key not available 
% see e.g.:
% https://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/
% must be run by PLDAPS_preprocessing
% 10-2020: CF updated for OEphys and newline

% 02-2021: similar to its predecessor, but searches specifically for .psort files
% and only downloads them and their associated continuous and event files

data = struct; % might as well start this here

isoephys = strcmp(remoteDir(end-5:end-1),'neuro');

% get file list from local dir
if ~exist(localDir,'file')
    mkdir(localDir);
end
allFiles = dir(localDir);
if length(allFiles)>2
    localFolderList = allFiles(3).name; % skip 1+2, they are are "." and ".."
    for n = 4:length(allFiles)
        localFolderList = [localFolderList newline allFiles(n).name];
    end
else
    localFolderList = [];
end

% get file list from remote dir
cmd = 'ifconfig en7 inet'; % check if on local MBI network or need VPN workaround
[~,ifstuff] = system(cmd);
if any(strfind(ifstuff,'172.')) % MBI machine
    IP = '172.30.3.33';
else  % probably off campus, try proxy IP (requires VPN)
    IP = '10.161.240.133'; 
end

cmd = ['ssh fetschlab@' IP ' find ' remoteDir ' -name "*.psort"'];
% [~,remoteFileList] = system(cmd, '-echo');
[~,remoteFileList] = system(cmd); % no echo

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
    % remoteFiles has the full path; grab the actual file name, and its containing folder
    slashes = strfind(remoteFiles{n},'/');
    remoteFileName = remoteFiles{n}(slashes(end)+1:end);
    remoteFolderName = remoteFiles{n}(slashes(end-1)+1:slashes(end)-1);
    remotePath = remoteFiles{n}(1:slashes(end));

    dateStart = strfind(remoteFileName,'20');
    if isempty(dateStart); continue; end
    thisDate = remoteFileName(dateStart(1):dateStart(1)+9); % dateStart(1) just in case '20' appears in a timestamp
    thisDate([5 8]) = []; % OpenEphys dates are always yyyy-mm-dd, so remove the dashes
    
    if any(strfind(dateStr,thisDate))

        remoteFileName
        
        if ~exist([localDir remoteFolderName],'file')
            mkdir([localDir remoteFolderName]);
        end
        if overwriteLocalFiles
            cmd1 = ['scp fetschlab@' IP ':' remoteFiles{n} ' ' localDir remoteFolderName '/'];
        else
            cmd1 = ['rsync -av --ignore-existing fetschlab@' IP ':' remoteFiles{n} ' ' localDir remoteFolderName '/'];
        end
        system(cmd1,'-echo');
        
        % find the original continuous file used for sorting
        currentFullfile = fullfile(localDir,remoteFolderName,remoteFileName);
        contfile_orig  = Psort_read_contfileOnly(currentFullfile);
        slash = strfind(contfile_orig,'\');
        if isempty(slash); slash = strfind(contfile_orig,'/'); end % (windows vs mac/linux)
        contFile = contfile_orig(slash(end)+1:end);

        dot = strfind(contFile,'.');
        underscore = strfind(contFile,'_');
        dot_underscore = [dot underscore];
        
        % channel num is found after CH, and before the next dot OR underscore
        startchan = strfind(contFile,'CH')+2;
        endchan = min(dot_underscore(dot_underscore>startchan))-1;
        chan = contFile(startchan:endchan);
        
        % we can use underscores again, to find the 'run' number
        if length(underscore)==1
            run = '1';
            evFile = 'messages.events';
        elseif length(underscore)==2
            run = contFile(underscore(2)+1:dot-1);
            evFile = ['messages_' run '.events'];
        else
            error('cannot parse contfile name: too many/few underscores');
        end
        
%         % or just grab all of them, in case if a mismatch
%         evFile = [remoteFiles{n}(1:slashes(end)) 'messages*.events'];

        if overwriteLocalFiles
            cmd2 = ['scp fetschlab@' IP ':' remotePath contFile ' ' localDir remoteFolderName '/'];
            cmd3 = ['scp fetschlab@' IP ':' remotePath evfile ' ' localDir remoteFolderName '/'];
        else
            cmd2 = ['rsync -av --ignore-existing fetschlab@' IP ':' remotePath contFile ' ' localDir remoteFolderName '/'];
            cmd3 = ['rsync -av --ignore-existing fetschlab@' IP ':' remotePath evFile ' ' localDir remoteFolderName '/'];
        end
        system(cmd2,'-echo');
        system(cmd3,'-echo');
        
        % save the filename, run, chan, and maybe other stuff
        data.remotefile{m,1} = remoteFiles{n};
        data.localfile{m,1} = currentFullfile;
        data.filename{m,1} = remoteFileName;
        data.chan(m,1) = str2double(chan);
        data.run(m,1) = str2double(run);
        
        m=m+1;
    end
end

if exist('currentFilename','var')
    % sort the data struct by date, and link up mapping and exp files
    keyboard
else
    disp('no matching files found');
end
disp('DONE.');