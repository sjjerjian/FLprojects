% createNexData

% SJ 06-14-2021
% process raw nexonar recorded text files into PLDAPS-style data struct 'nex'
% calls function processNexonarUDP on loaded .txt files
% afterwards, move BOTH resulting .mat file and .txt files to NAS

% 06-14-2021
% currently sending PLDAPs iTrial and unique rng for syncing with PLDAPS
% files in case of any mismatches...
% could also send unique_trial_number to be consistent with Trellis

% since nex can have issues with drops/repeats of trial info because of UDP
% packet drops, there is a subsequent cleanNexonar function which uses
% PLDAPS info anyway, so iTrial becomes useful

clear nex
 
subject = 'test';
dateRange = 20230216;
hdr.dataPath = ['C:\Users\fetschlab\data\' subject '\'];
% hdr.dataPath = ['C:\Users\fetschlab\data\'];

files = dir([hdr.dataPath '*.txt']);

dateStr = [];
for d = 1:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

for f = 1:length(files)
    
    if contains(files(f).name,subject)
        
        hdr.nexFilename = files(f).name;
        dateStart = strfind(hdr.nexFilename,'20');
        thisDate  = hdr.nexFilename(dateStart(1):dateStart(1)+7);
        
        if any(strfind(dateStr, thisDate))
            dot = strfind(hdr.nexFilename,'.');
            hdr.filename = [hdr.nexFilename(1:dot-1) '_nexonar.mat'];
        
%             data = load([hdr.dataPath hdr.nexFilename]);
            % SJ 09/2021, to handle variable column numbers
%             data = readtable([hdr.dataPath hdr.nexFilename],'ReadVariableNames',0,'HeaderLines',0,'Delimiter',' ');
%             data = table2array(data);

%           % SJ 01/2022
            data = textread([hdr.dataPath hdr.nexFilename]);
%             fid = fopen([hdr.dataPath hdr.nexFilename]);
%             data = textscan(fid,'%8.4f');
            
            if isempty(data)
                disp('no data found...skipping\n'),
                continue
            end
            
            fprintf('processing %s\n',hdr.nexFilename)
            nex = processNexonarUDP(data);
            
            save([hdr.dataPath hdr.filename],'nex','hdr')
        end
    end
end