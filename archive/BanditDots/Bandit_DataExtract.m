%Miguel Vivar-Lazo
%10/16/2018
%Fetsch Lab
%HK - adapted for the bandit data 1/22/18

%Folder = Path 'string'
%Subject = Name of subject 'string'
%dateRange = date of file 'number'
%setupfile = name of the file PLDAPS use for setting up the task 'string'
%time = time of file 'string' - HK removed
function [PDS2] = Bandit_DataExtract(folder, subject, dateRange, setupFile)
    % set up data structuremasterData
    PDS = struct;
 
    % search folder for matching files and extract the desired variables from PDS
    allFiles = dir(folder);
    for n = 3:length(allFiles) % skip 1+2, they are are "." and ".."
        for m = 1:length(dateRange)
            targetfile = [subject num2str(dateRange(m)) setupFile]; % HK - time removed
            if strfind(allFiles(n).name, targetfile) % if the target file is there,
                load([folder allFiles(n).name],'-mat'); % load it.
            end
        end
    end
    
    PDS2 = PDS; %Turn load file into struct variable
end

