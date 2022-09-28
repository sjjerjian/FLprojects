function au = dots3DMP_FRmatrix_fromDataStruct(dataStruct,par,timingInfo,conds,condlabels,opts)

% conds and condlabels together allow extraction of trial-averaged firing
% rates (within an interval and/or binned around a certain event), for any
% condition grouping

% TODO check that this works for choice/wager conditions

if ~isfield(opts,'keepMU'), opts.keepMU = 1; end
if ~isfield(opts,'calcTuning'), opts.calcTuning = 0; end
if ~isfield(opts,'smoothFR'), opts.smoothFR = 0; end
                
if opts.smoothFR==1 && ~isfield(opts,'convKernel')
    opts.convKernel = fspecial('average', [1 20]); % N bins wide boxcar (acausal, centered)
end

alignEvent  = timingInfo.alignEvent;
otherEvents = timingInfo.otherEvents;
tStart      = timingInfo.tStart;
tEnd        = timingInfo.tEnd;
binSize     = timingInfo.binSize;

%% now loop over the dataStruct

clear au
clear condlist
unitInd = 0;
for ses = 1:length(dataStruct) 
    if ~isfield(dataStruct(ses).data,par), continue, end
    temp = dataStruct(ses).data.(par);
    
    if opts.keepMU, unit_inds = find(temp.units.cluster_type<3);
    else,           unit_inds = find(temp.units.cluster_type==2); % SU only
    end

    Ntr = length(temp.events.trStart);
    condlist = nan(Ntr,length(condlabels));

    % get the list of conditions on each trial for *this* session
    for c=1:length(condlabels)
        if strcmp(par,'VesMapping') && strcmp(condlabels{c}(1:7),'heading')
            condlist(:,c) = cell2mat(temp.events.(condlabels{c}));
        else
            condlist(:,c) = temp.events.(condlabels{c});
        end
    end

    % (if we have heading as a condition), force very small ~eps headings to zero (some discrepancies)
    hdgCol = strcmp(condlabels,'heading');
    if sum(hdgCol)
        condlist(abs(condlist(:,hdgCol))<0.01,hdgCol) = 0;
    end

    % force ves to be lowest coh in the set, kluge
    cohCol = contains(condlabels,'coherence');
    modCol = strcmp(condlabels,'modality');
    if sum(cohCol) && sum(modCol)
        condlist(condlist(:,modCol)==1,cohCol) = min(condlist(:,cohCol));
    end

    % remove any unwanted trials - brfix, or condition not in the conds list
    isInCondList = false(1,Ntr);
    condI        = nan(1,Ntr); % SJ added 09/28/2022!!
    for itr=1:Ntr
        [isInCondList(itr),condI(itr)] = ismember(condlist(itr,:),conds,'rows');
    end
    goodtrs = temp.events.goodtrial & isInCondList;
    
    % condlist contains each unique condition (per row)
    % condI stores the condition of each trial by its index in conds
    condI    = condI(goodtrs);
    condlist = condlist(goodtrs,:);


    % startInd (+1) is the starting ind for units in the current session, unitInd will go from 1:length(unit_inds) for each alignment event
    startInd = unitInd;

    for iae=1:length(alignEvent)

        % for all of this, remember that alignment for time-res psth (and relative time of otherEvents) will
        % ALWAYS ALWAYS ALWAYS be to the first column in alignEvent

        % if just one entry is given, align to this event in events
        if length(alignEvent{iae})==1
            ae = temp.events.(alignEvent{iae}{1})';
        else % 2 events, then align to first one (tStart), and go up to tEnd from second one
             ae = [temp.events.(alignEvent{iae}{1}); temp.events.(alignEvent{iae}{2})]';
        end

        % if a third argument is given, re-assign ae as this proportion
        % away from the first event, towards the second event
        if length(alignEvent{iae})==3
            ae = ae(:,1) + (ae(:,2)-ae(:,1))*alignEvent{iae}{3};
        end
        ae = ae(goodtrs,:);

        % calculate times of 'otherEvents' relative to main alignment event
        oe_times = nan(size(ae,1),length(otherEvents{iae}));
        for ioe=1:length(otherEvents{iae})
            try
                oe_times(:,ioe) = temp.events.(otherEvents{iae}{ioe})(goodtrs)' - ae(:,1);
            catch
%                 fprintf('Could not calculate relative time of %s
%                 event\n',otherEvents{iae}{ioe});
            end
        end

        unitInd = startInd; % back to the 'top' (i.e. startInd), so that we are adding to the same unit index for each alignment event
        
        for u = 1:length(unit_inds)
            unitInd = unitInd+1; % ADD 1 to the INDEX HERE

            au.hdr.unitDate(unitInd) = dataStruct(ses).date;
            au.hdr.unitSet(unitInd)  = dataStruct(ses).set;
            au.hdr.unitID(unitInd)   = temp.units.cluster_id(unit_inds(u));
            au.hdr.unitType(unitInd) = temp.units.cluster_type(unit_inds(u));


            % calculate time-resolved and average firing rate within desired intervals
            [fr, x, fr_mean, durs, aeI]         = trial_psth(temp.units.spiketimes{unit_inds(u)},ae,'tStart',tStart(iae),'tEnd',tEnd(iae),'binSize',binSize);

            % now average activity across trials of each condition in conds

%             if unitInd == 352, keyboard, end
            condFR      = nan(size(conds,1),size(fr,2));
            condFR_mean = nan(size(conds,1),size(fr_mean,2));
            condFR_sem  = condFR_mean;
            durs_median = condFR_mean;

            for ic = 1:size(conds,1)
                
                condFR(ic,:)     = nanmean(fr(condI==ic,:),1);
%                 condFR_sem(ic,:) = nanstd(fr(condI==ic,:),[],1)/sqrt(sum(condI==ic));

                condFR_mean(ic)  = nanmean(fr_mean(condI==ic));
                condFR_sem(ic)   = nanstd(fr_mean(condI==ic))/sqrt(sum(condI==ic));

                % store median time of 'otherEvents' relative to alignment event, for each condition (for marking on PSTHs)
                au.times.evTimes{iae}(ic,:,unitInd) = nanmedian(oe_times(condI==ic,:));

                durs_median(ic)  = nanmedian(durs(condI==c)); % and median duration

            end


            if opts.calcTuning && sum(hdgCol)

                % pull out headings list and trial headings
                hdgs = unique(conds(:,hdgCol));
                hdg  = condlist(:,hdgCol);

                [uconds,~,ic2] = unique(conds(:,~hdgCol),'rows');
                for uc=1:size(uconds,1)
                    theseTrials = all(condlist(:,~hdgCol)==uconds(uc,:),2);
                    thisCond_trialFR = fr_mean(theseTrials)';
                    thisCond_meanFR  = condFR_mean(ic2==uc);

                    au.stat.ntrsCond(uc,unitInd) = sum(theseTrials);

                    % should keep track of the strength of pref too
                    % (sum(Right)-sum(Left)/sum(Right)+sum(Left)
                    % empirical heading pref - which heading dir (left or right) has higher firing?
                    if sum(theseTrials)>0
                        au.stat.prefAmp(uc,unitInd) = (sum(thisCond_meanFR(hdgs>0)) - sum(thisCond_meanFR(hdgs<0))) / sum(thisCond_meanFR(hdgs~=0));
                        au.stat.prefDir(uc,unitInd) = double(sum(thisCond_meanFR(hdgs>0)) > sum(thisCond_meanFR(hdgs<0)))+1; % +1 since choice is 1,2 for L,R
                    else
                        au.stat.prefAmp(uc,unitInd) = NaN;
                        au.stat.prefDir(uc,unitInd) = NaN;
                    end

                    % or where is the peak?
                    %                 [~,peak] = max(thisCond_meanFR);
                    %                 au.prefDir(uc,unitInd) = sign(hdgs(peak));

                    % or fit von Mises (work in progress)
                    % vonMises params: b = [ampl, kappa, theta, offset]
                    %                 guess = [max(thisCond_meanFR)-min(thisCond_meanFR), 1.5, peak, min(thisCond_meanFR)];
                    %                 [beta, fval] = fminsearch(@(x) tuning_vonMises_err(x,hdg(theseTrials),thisCond_trialFR), guess);

                    % simple calc of tuning significance
                    [au.stat.pHdg(uc,unitInd)]=anova1(thisCond_trialFR,hdg(theseTrials),'off');

                end
            end
            
            if opts.smoothFR
                condFR = smoothRaster ( condFR , opts.convKernel );
            end

            % these could be slightly different lengths across units
            % so keep as cells and then we will 'matricise' them outside
            % the loop. maybe there is a more efficient way to do this
            temp_xvec{iae}{unitInd} = x;
            temp_psth{iae}{unitInd} = condFR;
            
            au.data.muFRs{iae}(:,unitInd) = condFR_mean;
            au.data.seFRs{iae}(:,unitInd) = condFR_sem;
            au.data.durs{iae}(:,unitInd)  = durs_median;       
        end
    end
end

% finally, transfer the firing rate of each unit for each condition into a
% big matrix (conds x timepoints x units)
% set the start and end points using whichever xvec is the longest
% ... in practice, there might just be one or two bins difference on either
% side across sessions, and this is only really necessary if alignEvent
% involves two events with variable interval between them

for iae=1:length(alignEvent)
    [maxlenx,pos] = max(cellfun(@length, temp_xvec{iae}));

    au.data.PSTHs{iae} = nan(size(temp_psth{iae}{1},1),maxlenx,unitInd);
    au.times.xvec{iae} = temp_xvec{iae}{pos};

    for u=1:unitInd
        this_psth = temp_psth{iae}{u};
        if aeI == 1
            st = 1;
            en = size(this_psth,2);
        elseif aeI == 2
            st = maxlenx - size(this_psth,2) + 1;
            en = maxlenx;
        end
        
        au.data.PSTHs{iae}(:,st:en,u) = this_psth;

    end
end
    

