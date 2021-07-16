function nsdata = createDataStructure_NS(subject,date,filenum,pldaps_filetime,filepath)
% requires NeuroShare toolbox and custom getdata_NS function

% creates data structure comparable to PLDAPS, with trial info, outcomes,
% and key event times - one entry per trial


% addpath(genpath('C:\Program Files (x86)\Ripple\Trellis\Tools'));
% addpath(genpath('C:\Users\fetschlab\Documents\MATLAB\dots3DMP'));
                    
% these should be the function arguments
% subject  = 'lucio';
% date     = 20210517;      
% filenum  = 1;             % filename in Trellis (auto-increment numbers)
% pldaps_filetime = 1340;   % PLDAPS file names have start time at end 
% filepath = 'C:\Users\fetschlab\Trellis\dataFiles\';

pldaps_filename = sprintf('%s%ddots3DMP%d',subject,date,pldaps_filetime);
filename = sprintf('%s%ddots3DMP%0.4d.nev',subject,date,filenum);
completeFilePath = fullfile(filepath,filename);

nsdata.pldaps_filename = pldaps_filename;
nsdata.trellis_filename = filename;
nsdata.trellis_filepath = completeFilePath;

dataType    = 'Digital';
dataChannel = 'parallel';  
[nsEvents] = getdata_NS(completeFilePath, dataType, dataChannel);
tr_start = strfind(nsEvents.data>=2^9,[1 1 1 1 1 1 1]); % find clocktime for trial delimiter

infolabels = {'headingI','modality','coherenceI','deltaI','choice','correct','PDW'};

%%
for t=1:length(tr_start)
    
    if t==length(tr_start)
        t_en = length(nsEvents.data);
    else
        t_en = tr_start(t+1)-1;
    end
    
    tr_events = nsEvents.data(tr_start(t):t_en);
    tr_times  = nsEvents.time(tr_start(t):t_en);
    
    tr_info = tr_events(tr_events>=2^8 & tr_events<2^9) - 2^8;
   
    nsdata.pldaps_iTr(t) = tr_events(7)-(2^9+2^8);
    nsdata.pldaps_uTr{t} = tr_events(1:6)-(2^8);
    
    
    % try to be sure that each tr_events only contains ONE trial's events,
    % i.e. there should only be two 1 strobes for start and end of PLDAPS
    % trial.
    % PLDAPS trial numbers should also be incrementing by 1
    if (t>1 && nsdata.pldaps_iTr(t)-nsdata.pldaps_iTr(t-1)~=1) || sum(tr_events==1)~=2
        keyboard
    end
    
    for i=1:7%length(tr_info)
        try
            nsdata.(infolabels{i})(t)  = tr_info(i);
        catch
            nsdata.(infolabels{i})(t)  = nan;
        end
    end
   
    try nsdata.trStart(t)   	= tr_times(find(tr_events==1,1)); 
    catch, nsdata.trStart(t)    = nan; 
    end
    try nsdata.trEnd(t)         = tr_times(find(tr_events==1,1,'last')); 
    catch, nsdata.trEnd(t)      = nan; end
    try nsdata.fpOn(t)          = tr_times(tr_events==2); 
    catch, nsdata.fpOn(t)       = nan; end
    try nsdata.fix(t)           = tr_times(tr_events==3); 
    catch, nsdata.fix(t)        = nan; end
    try nsdata.targsOn(t)       = tr_times(tr_events==4); 
    catch, nsdata.targsOn(t)    = nan; end
    try nsdata.stimOn(t)        = tr_times(find(tr_events==5,1)); 
    catch, nsdata.stimOn(t)     = nan; end
    try nsdata.stimOff(t)       = tr_times(find(tr_events==5,1,'last')); 
    catch, nsdata.stimOff(t)    = nan; end
    try nsdata.sacc(t)          = tr_times(tr_events==6); % 'RT'
    catch, nsdata.sacc(t)       = nan; end
    try nsdata.targHold(t)    	= tr_times(tr_events==7); % CHOICE
    catch, nsdata.targHold(t)   = nan; end
    try nsdata.posttargHold(t)  = tr_times(tr_events==8); % PDW
    catch, nsdata.posttargHold(t) = nan; end
    try nsdata.reward(t)        = tr_times(tr_events==9); 
    catch, nsdata.reward(t)     = nan; end
end


