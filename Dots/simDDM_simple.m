%% simDDM_simple.m
% simple simulation of a 1D drift-diffusion model in the context of a
% motion discrimination task (see e.g. Shadlen et al. 2006)
% CF circa 2015

clear all; close all;

ntrials = 5000;

% different levels of motion strength ('coherence')
cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 -eps eps 0.032 0.064 0.128 0.256 0.512];

% delta-T (time step, in ms)
dT = 1;
% max duration of stimulus
maxdur = 2000;
timeAxis = 0:dT:maxdur;

% randomize coherence and duration for all trials (n draws with replacement)
coh = randsample(cohs,ntrials,'true')';



%% model parameters (play around with these to see their effects)
k = 0.3; % 'drift rate' or sensitivity term: a constant converting stimulus
         % strength into units of momentary evidence
sigma = 1; % standard deviation of momentary evidence; often fixed at 1
B = 25; % height of the bound, or threshold, for decision termination
Tnd = 300; %non-decision time

%% simulate the diffusion process

% initialize variables
dv = nan(ntrials,length(timeAxis)); % the DV as a function of n and time
choice = nan(ntrials,1); % vector to store choices: 2=right=positive, 1=left=negative
RT = nan(ntrials,1); % reaction time
finalV = nan(ntrials,1); % endpoint of the DV
hitBound = nan(ntrials,1);

tic
for n = 1:ntrials
    % mean of momentary evidence
    mu = k * coh(n);

    % diffusion process: slow version***
    if n==1 % ***run this once, for illustration purposes
        momentaryEvidence = nan(maxdur,1);
        for t = 1:maxdur
            momentaryEvidence(t) = randn*sigma + mu; % random sample with mean=mu, s.d.=sigma
        end
        dv(n,1) = 0;
        dv(n,2:maxdur+1) = cumsum(momentaryEvidence);
        figure; plot(dv(n,:)); hold on; 
        plot(1:length(dv),ones(1,length(dv))*B,'g-');
        plot(1:length(dv),ones(1,length(dv))*-B,'r-');
        xlabel('Time (ms)'); ylabel('Accum. evidence (DV)');
        % (evidence is shown continuing to accumulate past the bound,
        % although it realy stops there; this can be useful
        % for diagnosing problems with certain parameter settings)
    end    

    % faster version: does not require a FOR loop over the variable t
    dv(n,:) = [0, cumsum(normrnd(mu,sigma,1,maxdur))];
        
    tempRT = find(abs(dv(n,:))>=B, 1);
    if isempty(tempRT) % did not hit bound
        RT(n) = maxdur;
        finalV(n) = dv(n,RT(n));
        hitBound(n) = 0;
    else % hit bound
        RT(n) = tempRT;
        finalV(n) = B*sign(dv(n,RT(n)));
        hitBound(n) = 1;
    end
    choice(n) = sign(finalV(n));  
end
toc

RT = RT + Tnd;

% quick sanity check to see if params give reasonable performance
pCorrect_total = sum(sign(choice)==sign(coh)) / ntrials


%% format data as in experimental data files and generate output structs

coh(coh==0) = sign(randn)*eps; % should have no actual zeros, but if so, sign them randomly;
                               % this is just to assign a direction and correct/error
data.correct = choice==sign(coh);
data.direction = nan(ntrials,1);
data.direction(coh>0) = 0;
data.direction(coh<0) = 180;
coh(abs(coh)<1e-6) = 0; % now go back to one 'zero'
data.coherence = abs(coh);
data.scoh = coh;

data.choice = choice;
data.choice(data.choice==-1) = 0; % code elsewhere assumes 0s and 1s
data.RT = RT/1000; % convert to seconds

conftask = 0; RTtask = 1; RTCorrOnly = 0;
parsedData = Dots_parseData(data,conftask,RTtask,RTCorrOnly);

% plot
cohs = unique(coh); wFit = 0; forTalk = 0;
Dots_plot(parsedData,cohs,conftask,RTtask,wFit,forTalk)



