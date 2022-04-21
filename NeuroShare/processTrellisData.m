function processTrellisData(info)

% process Trellis Data
%
% extract behavioural events, and online detected/sorted spikes
% extract and save raw neural recordings in binary format for kilosort
% sorting

% SJ updated 01/2022

% This is the main script to run after a recording for extracting and
% formatting synced task and neural data!

% 1.
% load .nev file recorded via Trellis, and extract task stimulus and
% behavioral events/times in 'Trellis time', for aligning neural data to.
% To store header information regarding matching PLDAPS file, manually
% enter list of pldaps_filetimes below - list should correspond exactly to
% number of .nev files. Stem of pldaps_filename should have been set already in
% Trellis recording save field

% 2.
% for the same nevFiles, also extract  online threshold crossing and
% spikes, if desired.

% 3. extract raw neural signals, create .bin files for Kilosort

% REQUIRES
% - neuroshare toolbox
%       - getdata_NS to extract raw data
%       - createEventStruct_NS to re-organize raw data events into trials
% - DBCchan_remap to get correct channels for extracting spikes

% clear;clc

%%

% files should be in date directory
nevFiles = dir([info.filepath '*.nev']);

if isempty(nevFiles)
    disp('No files found...check folder')
end

for f=1:length(nevFiles)
%     if f==2,keyboard,end

    % ASSUME THAT EVENT FILES/PLDAPS FILES ARE MATCHED FOR SIMULTANEOUS
    % RECORDINGS
    if contains(nevFiles(f).name,info.subject) && ~isnan(info.pldaps_filetimes(f))

        %         dateStart = strfind(nevFiles(f).name,'20');
        %         thisDate  = nevFiles(f).name(dateStart(1):dateStart(1)+7);

        dot = strfind(nevFiles(f).name,'.nev');
%         filenum = str2double(nevFiles(f).name(dot-4:dot-1));

        pldaps_filename = sprintf('%s%d%s%04d',info.subject,info.date,info.par{f},info.pldaps_filetimes(f));
        bin_filename = nevFiles(f).name(1:dot-1);
        savename = [bin_filename '_RippleEvents.mat'];

        fprintf('\n processing task events for file %s [%s] (%d of %d)...\n',nevFiles(f).name,sprintf('%04d',info.pldaps_filetimes(f)),f,length(nevFiles))
        nsEvents = createEventStruct_dots3DMP_NS(nevFiles(f).name,info.par{f},pldaps_filename,info.filepath); % newer 02/2022
        %         nsEvents = createEventStruct_NS(subject,thisDate,paradigm,filenum,pldaps_filetime,filepath); 

        allRawData = [];
        for pb = 1:length(info.probe_type) % loop over probes/electrodes

            if ~isempty(info.chanlist{pb})
%                 completeFilePath = fullfile(info.filepath,nevFiles(f).name);
                for ch=1:length(info.chanlist{pb})
                    thisChan = info.chanlist{pb}(ch);
 
                    % worth getting online sorted data for chanInterest?
                    if ismember(info.chanlist{pb}(ch),info.chanInterest{pb}{f})
                        fprintf('extracting online spikes from ch%d\n',thisChan)
                        nsEvents.spkData.data{ch} = getdata_NS(nsEvents.hdr.Ripple_filepath, 'Spike', thisChan);
                        nsEvents.spkData.chs(ch) = thisChan; 
                    end

                    % get raw continuous data from all channels
                    %             fprintf('processing raw data from chan %d\n',thisChan);
                    [~,name,~] = fileparts(nevFiles(f).name);
                    temp = getdata_NS(fullfile(info.filepath,[name '.ns5']), 'Raw', thisChan);
                    rawData(:,ch) = temp.analogData';
                end
            end
            % append each probe to a full data matrix
            allRawData = cat(2,rawData,allRawData);

            % write rawData to int16 binary files for Kilosort
            fid = fopen([info.filepath bin_filename '.bin'],'w');
            fwrite(fid,int16(allRawData),'int16');
            fclose(fid);

            % for debugging
            %             fid = fopen([filepath bin_filename '_ch' num2str(thisChan) '.bin'],'r');
            %             X = fread(fid,'int16');
            %             fclose(fid);

            % save Events
            save([info.filepath savename],'nsEvents');

        end
    end
end








