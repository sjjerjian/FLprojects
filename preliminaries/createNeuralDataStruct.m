%% createNeuralDataStruct;

% grabs some neural data from server and creates a dataCell containing

% requires password(s) to be entered if a auth key not available
% see e.g.:
% https://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/

% 03-2022 SJ, separate to PLDAPS preprocessing, for spike sorted waveforms,
% timestamps, and RippleEvents

cmd = 'ifconfig en7 inet'; % check if on local MBI network or need VPN workaround
[~,ifstuff] = system(cmd);
if ~useVPN || any(strfind(ifstuff,'172.'))
    IPadd = 'fetschlab@172.30.3.33'; % MBI machine
else
    IPadd = 'fetschlab@10.161.240.133'; % probably off campus, try proxy IP (requires VPN)
end

% get file list from local dir
if ~exist(localDir,'dir')
    mkdir(localDir);
end
allFiles = dir(localDir);
if length(allFiles)>2
    localFileList = allFiles(3).name; % skip 1+2, they are are "." and ".."
    for n = 4:length(allFiles)
        localFileList = [localFileList newline allFiles(n).name];
    end
else
    localFileList = '';
end

% get folder list of recording sessions from remote dir

cmd = ['ssh ' IPadd ' ls ' remoteDir];
[~,remoteFolderList] = system(cmd, '-echo');
% [~,remoteFolderList] = system(cmd);
if any(strfind(remoteFolderList,'timed out')); error(remoteFolderList); end


% check each folder for a match to the desired date range, and pull in the
% RippleEvents and info files from the server

newlines = strfind(remoteFolderList,newline);
remoteFolders = cell(length(newlines),1);
m = 1; m2 = 1;
sess = 0;
for n = 1:length(newlines)
    if n==1
        remoteFolders{n} = remoteFolderList(1:newlines(1)-1);
    else
        remoteFolders{n} = remoteFolderList(newlines(n-1)+1:newlines(n)-1);
    end
    dateStart = strfind(remoteFolders{n},'20');
    if isempty(dateStart) || contains(remoteFolders{n},'Impedance'); continue; end
    
    thisDate = remoteFolders{n}(dateStart(1):dateStart(1)+7); % dateStart(1) just in case '20' appears in a timestamp
    
    if any(strfind(dateStr,thisDate))
        currentFolderList{m,1} = remoteFolders{n}; m=m+1;
        
        cmd = ['ssh ' IPadd ' ls ' remoteDir remoteFolders{n}];
        %         [~,thisFileList] = system(cmd, '-echo');
        [~,thisFileList] = system(cmd);
        
        newlines2 = strfind(thisFileList,newline);
        
        matFiles = cell(length(newlines2),1);
        for n2 = 1:length(newlines2)
            if n2==1
                matFiles{n2} = thisFileList(1:newlines2(1)-1);
            else
                matFiles{n2} = thisFileList(newlines2(n2-1)+1:newlines2(n2)-1);
            end
            [~,~,ext] = fileparts(matFiles{n2});
            if ~strcmp(ext,'.mat');  continue; end
            currentFileList{m2,1} = matFiles{n2}; m2=m2+1;
            
            if ~contains(localFileList,matFiles{n2}) || overwriteLocalFiles % always copy if overwrite option selected
                cmd = ['scp -r ' IPadd ':' remoteDir remoteFolders{n} '/' matFiles{n2} ' ' localDir];
                system(cmd,'-echo');
            else
                disp([remoteFolders{n} '/' matFiles{n2} ' exists locally, not copied']);
            end
            
        end
    end
end

%% create dataCell

% maybe separate this and load in existing datacell to append to?
sess = 1;
m=1;
clear currentFileList dataCell
dataCell = struct();
for n = 1:length(currentFolderList)
    disp(currentFolderList{n})
    if isempty(strfind(currentFolderList{n},'20')) || contains(currentFolderList{n},'Impedance'); continue; end
    
    clear info
    load(fullfile(localDir,[subject currentFolderList{n} 'dots3DMP_info.mat']));
    
    
    for f=1:length(info.trellis_filenums)
        
        if strcmp(info.par{f},paradigm) && ~isnan(info.pldaps_filetimes(f))
            
            PDSfilename = [info.subject num2str(info.date) info.par{f} num2str(info.pldaps_filetimes(f))];
            NSfilename  = sprintf('%s%ddots3DMP%04d_RippleEvents.mat',info.subject,info.date,info.trellis_filenums(f));
            
            fprintf('adding units from %s\n',NSfilename)
            load(fullfile(PDSdir,PDSfilename));
            load(fullfile(localDir,NSfilename));
            
            % cross-ref with PDS struct, insert the heading/coh/delta inds
            nsEvents = nsEventConditions(nsEvents,PDS);
            
            dataCell(sess).nsEvents    = nsEvents;
            dataCell(sess).paradigm    = info.par{f};
            dataCell(sess).date        = info.date;
            dataCell(sess).fileID      = info.trellis_filenums(f);
            
            dataCell(sess).PDSfilename = PDSfilename;
            dataCell(sess).NSfilename  = NSfilename;
            dataCell(sess).rec_group   = info.rec_group;
            
            % need to access the sorted data - it's in a subfolder of local dir
            % if info.probe_type 'single electrode'
            if contains(info.probe_type{1},'Single')
                remoteDirSpikes = sprintf('/var/services/homes/fetschlab/data/%s/%s_neuro/%d/%s%ddots3DMP%04d/',subject,subject,info.date,subject,info.date,info.trellis_filenums(f));
                mountDir = sprintf('/Volumes/homes/fetschlab/data/%s/%s_neuro/%d/%s%ddots3DMP%04d/',subject,subject,info.date,subject,info.date,info.trellis_filenums(f));
            else
                remoteDirSpikes = sprintf('/var/services/homes/fetschlab/data/%s/%s_neuro/%d/%s%d_%d/',subject,subject,info.date,subject,info.date,info.rec_group(f));
                mountDir = sprintf('/Volumes/homes/fetschlab/data/%s/%s_neuro/%d/%s%d_%d/',subject,subject,info.date,subject,info.date,info.rec_group(f));
            end
            
            cmd = ['ssh ' IPadd ' ls ' remoteDirSpikes];
            [~,remoteFileList] = system(cmd);
            
            newlines = strfind(remoteFileList,newline);
            remoteFiles = cell(length(newlines),1);
            
            for nl = 1:length(newlines)
                if nl==1
                    remoteFiles{nl} = remoteFileList(1:newlines(1)-1);
                else
                    remoteFiles{nl} = remoteFileList(newlines(nl-1)+1:newlines(nl)-1);
                end
                
                % MKSORT (single electrode recordings)
                if contains(remoteFiles{nl},'waveforms')
                    currentFileList{m,1} = remoteFiles{nl}; m=m+1;
                    load([mountDir remoteFiles{nl}]);
                    
                    uunits = unique(waveforms.units);
                    
                    for u=1:length(uunits)
                        if uunits(u)>0
                            dataCell(sess).spikeTimes{u-1} = waveforms.spikeTimes(waveforms.units==uunits(u));
                            dataCell(sess).waves{u-1}      = waveforms.alignedWaves(:,waveforms.units==uunits(u));
                            dataCell(sess).unitID(u-1)     = uunits(u);
                            %                             dataCell(sess).unitRating(u-1) = waveforms.ratings.ratings(uunits(u));
                            dataCell(sess).numSpikes(u-1)  = sum(waveforms.units==uunits(u));
                            
                        end
                    end
                    
                    % KILOSORT
                elseif contains(remoteFiles{nl},'amplitudes.npy')
                    currentFileList{m,1} = remoteFiles{nl}; m=m+1;
                    sp = loadKSdir(mountDir);
                    
                    
                    % shift the spike times from Kilosort because of the
                    % concatenated recordings, back to 'time 0' for the
                    % individual nsEvents, so that we can align to task
                    % events
                    timeStampsShifted = nsEvents.analogInfo.timeStampsShifted ./ double(nsEvents.analogInfo.Fs);
                    thisBlockSpikes   = sp.st >= timeStampsShifted(1) & sp.st < timeStampsShifted(2);
                    
                    % NOTE: spike templates is only equal to spike clusters
                    % before any manual curation, then spike clusters will
                    % change to reflect new assignments, spike templates is
                    % fixed. so used spike clusters for assigments here!
                    
                    for u=1:length(sp.cids)
                        if sp.cgs(u)<3 % 3 denotes noise clusters, skip
                            theseSpikes = sp.clu==sp.cids(u) & thisBlockSpikes;
                            dataCell(sess).spikeTimes{u}    = sp.st(theseSpikes) - timeStampsShifted(1);
                            dataCell(sess).clusterID(u)     = sp.cids(u);
                            dataCell(sess).clusterGroup(u)  = sp.cgs(u); % 1 - MU, % 2 - SU %
                            dataCell(sess).numSpikes(u)     = sum(theseSpikes);
                            % should get the waveforms out too? other
                            % summary info?
                        end
                    end
                end
            end
        end
        sess=sess+1;
    end
end

disp('done.');
