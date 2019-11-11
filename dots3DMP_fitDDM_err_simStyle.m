function [err,data] = dots3DMP_fitDDM_err_simStyle(param, data, options)

global call_num

keyboard

%     param = getParam(param,guess,fixed); % ignore this for now
kves = param(1);
kvis = param(2);
B = abs(param(3)); % don't accept negative bound heights
sigma = 1;

ntrials = length(data.heading);

cohs = unique(data.coherence); % visual coherence levels
hdgs = [-10 -5 -2.5 -1.25 -eps eps 1.25 2.5 5 10]; % heading angles
                            % (map fn seems to require even number of diff 
                            % levels (or just no zero), so we use +/- eps)
deltas = unique(data.delta); % conflict angle; positive means vis to the right
mods = [1 2 3]; % stimulus modalities: ves, vis, comb
duration = 2000; % stimulus duration (ms)

theta = 0.6; % in old models, threshold for sure-bet choice in units of log odds correct.
             % not used for analog confidence judgment, but needed for the makeLogOddsCorrMap func

maxdur = duration;
k = kves*(1/3) + (kvis*cohs(1))*(1/3) + (kvis*cohs(2))*(1/3);

[logOddsMapR, logOddsMapL, logOddsCorrMap, tAxis, vAxis] = makeLogOddsCorrMap_3DMP(hdgs,k,B,sigma,theta,maxdur,1);
% uses Fokker-Planck equation to propagate the probability density of the DV,
% as in Kiani & Shadlen 2009. Required for readout of confidence, although
% a simpler heuristic could be used (conf proportional to accum evidence)

dur = ones(ntrials,1) * duration;


%% bounded evidence accumulation

% assume momentary evidence is proportional to sin(heading),
% as in drugowitsch et al 2014

dv = nan(ntrials,max(dur));
choice = nan(ntrials,1);
RT = nan(ntrials,1);
finalV = nan(ntrials,1);
hitBound = zeros(1,ntrials);
logOddsCorr = nan(ntrials,1);
expectedPctCorr = nan(ntrials,1);
conf = nan(ntrials,1);

hdg = data.heading;
coh = data.coherence;
delta = data.delta;
modality = data.modality;

tic
for n = 1:ntrials

%     n % progress check, can remove
    
    % slow version
    x = nan(1,dur(n));
    for t = 1:dur(n)
        switch modality(n)
            case 1 % ves
                mu = kves * sind(hdg(n)); % mean of momentary evidence
                x(t) = mu + randn*sigma; % random sample with mean=mu, s.d.=sigma
            case 2 % vis
                mu = kvis*coh(n) * sind(hdg(n));
                x(t) = mu + randn*sigma;
            case 3 % comb
                % positive delta defined as ves to the left, vis to the right
                muVes = kves        * sind(hdg(n) - delta(n)/2);
                muVis = kvis*coh(n) * sind(hdg(n) + delta(n)/2);

                wVis = sqrt( (kvis*coh(n))^2 / ((kvis*coh(n))^2 + kves^2) );
                wVes = sqrt( kves^2          / ((kvis*coh(n))^2 + kves^2) );

                x(t) = wVis*(muVis+randn*sigma) + wVes*(muVes+randn*sigma);
        end
    end
    dv(n,1) = 0; % initial condition
    dv(n,2:dur(n)+1) = cumsum(x);
    if n==1 % plot one trial as sanity check
        figure(10); plot(dv(n,:)); hold on; 
        plot(1:size(dv,2),ones(1,size(dv,2))*B,'g-');
        plot(1:size(dv,2),ones(1,size(dv,2))*-B,'r-');
        xlabel('Time (ms)'); ylabel('Accum. evidence (DV)');
        ylim([-B-5 B+5]);
        % (evidence is shown continuing to accumulate
        % past the bound, although it realy stops there)
    end
    

%     % faster version: avoids FOR loop over the variable t
                    %     [work in progress]
%     switch stimtype(n)
%         case 1
% %             x(t) = muVes + randn*sigma;
% 
%             dv(n,:) = [0, cumsum(normrnd(muVes,sigma,1,dur(n)))];
%         case 2
% %             x(t) = muVis + randn*sigma;
%             dv(n,:) = [0, cumsum(normrnd(muVis,sigma,1,dur(n)))];
%         case 3
% %             x(t) = wVis*(muVis+randn*sigma) + wVes*(muVes+randn*sigma);
% %             mu = wVis*(muVis+randn*sigma) + wVes*(muVes+randn*sigma);
%             % OR
%             mu = wVis*muVis + wVes*muVes; % (and may need to increase sigma)
%             dv(n,:) = [0, cumsum(normrnd(mu,sigma,1,dur(n)))];
%     end
        

    cRT = find(abs(dv(n,:))>=B, 1);
    if isempty(cRT) % did not hit bound
        RT(n) = max(dur);
        finalV(n) = dv(n,max(dur));
        hitBound(n) = 0;
    else % hit bound
        RT(n) = cRT;
        finalV(n) = B*sign(dv(n,cRT));
        hitBound(n) = 1;
    end    
    choice(n) = sign(finalV(n));
    
    % use map to look up log-odds that the motion is rightward
    diffV = abs(vAxis-finalV(n));
    diffT = abs(tAxis-RT(n));
        
    thisV = find(diffV==min(diffV));
    thisT = find(diffT==min(diffT));
    logOddsCorr(n) = logOddsCorrMap(thisV(1), thisT(1));
    expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
    conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1

end
toc

choice(choice==1) = 2; choice(choice==-1) = 1; % 1=left, 2=right

% calculate error


% RESUME







end


