
% SJ 11/2021
% testing

% compare time delays between events in PDS struct and events in NS...
% seems to be minimal (thanks digital strobes!)
% TODO:
% reassign stimulus condition index lists according to PDS conditions
% 
% dealing with spikes data, woohoo!

ntrs = length(PDS.data);

trstarts_PLDAPS = nan(ntrs,1);
trstarts_Ripple = nan(ntrs,1);


for t=1:ntrs
    
    data = PDS.data{t};
    stimConditions = PDS.conditions{t}.stimulus;
    
%     trstarts_PLDAPS(t,1) = data.trstart;
    trstarts_PLDAPS(t) = data.datapixx.unique_trial_time(1);

    trstarts_Ripple(t,1) = nsEvents.EventTimes.trStart(t);
end

trstarts_PLDAPS = trstarts_PLDAPS - trstarts_PLDAPS(1) + trstarts_Ripple(1);
fo = fit(trstarts_Ripple,trstarts_PLDAPS,'poly1');
trstarts_fixed = fo(trstarts_Ripple);

% SJ 11/2021
% <1ms lag between Ripple event times and PLDAPS times under assumed linear fit.
% can correct if we want by shifting according to fit...

% fnames = fieldnames(nsEvents.EventTimes);
% for f=1:length(fnames)
%     nsEvents.EventTimes.(fnames{f}) = fo(nsEvents.EventTimes.(fnames{f}));
% end



 