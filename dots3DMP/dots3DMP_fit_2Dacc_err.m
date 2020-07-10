function [err,fit] = dots3DMP_fit_2Dacc_err(param, data, options)

global call_num

ntrials = length(data.heading);

mods = unique(data.modality);
cohs = unique(data.coherence); % visual coherence levels
hdgs = unique(data.heading);
deltas = unique(data.delta);

duration = 2000; % stimulus duration (ms)
dur = ones(ntrials,1) * duration;

ks    = param(1);
sigma = param(2);
B     = abs(param(3)); % don't accept negative bound heights

kves = ks;
kvis = [.666 1.5]*ks;

sigmaVes = sigma; % std of momentary evidence
sigmaVis = [sigma sigma]; % allow for separate sigmas for condition, coherence

Tnd = 0.4; % non-decision time (or make this the mean of a dist?)

maxdur = duration;
% assume the mapping is based on an equal amount of experience with the 
% *three* levels of reliability (ves, vis-low, vis-high) hence k and sigma
% are their averages
k = mean([kves kvis]);

% [~, ~, logOddsCorrMap, tAxis, vAxis] = makeLogOddsCorrMap_3DMP(hdgs,k,B,sigma,maxdur,0);
% uses Fokker-Planck equation to propagate the probability density of the DV,
% as in Kiani & Shadlen 2009. Required for readout of confidence, although
% a simpler, non time-dependent relationship (conf proportional to accum
% evidence) could be used

R.t = 0.001:0.001:duration/1000;
R.Bup = B;
R.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
R.lose_flag = 1;
R.plotflag = 0; % 1 = plot, 2 = plot and export_fig

P =  images_dtb_2d(R);

% create acceleration and velocity profiles (arbitrary for now)
% SJ 04/2020
% Hou et al. 2019, peak vel = 0.37m/s, SD = 210ms
vel = normpdf(1:duration,duration/2,210);
vel = 0.37*vel./max(vel);
acc = gradient(vel)*1000; % multiply by 1000 to get from m/s/ms to m/s/s

% normalize
vel = vel./max(vel);
acc = abs(acc./max(acc)); % (and abs)


%% bounded evidence accumulation

choice = nan(ntrials,1);
RT = nan(ntrials,1);
finalV = nan(ntrials,1);
hitBound = zeros(1,ntrials);
logOddsCorr = nan(ntrials,1);
expectedPctCorr = nan(ntrials,1);
conf = nan(ntrials,1);

modality = data.modality;
hdg = data.heading;
coh = data.coherence;
delta = data.delta;

% ME is now a draw from bivariate normal with mean vector Mu and covariance matrix V
% start with correlation matrix:
S = [1 -1/sqrt(2) ; -1/sqrt(2) 1];
% -1/sqrt(2) is the correlation for our version of images_dtb;
% can only be changed with an update to the flux file

for n = 1:ntrials

    switch modality(n)
        case 1
            mu = kves * sind(hdg(n)) / 1000; % mean of momentary evidence
                % (I'm guessing drift rate in images_dtb is per second, hence div by 1000)
            s = [sigmaVes sigmaVes]; % standard deviaton vector (see below)
        case 2
            mu = kvis(cohs==coh(n)) * sind(hdg(n)) / 1000;
            s = [sigmaVis(cohs==coh(n)) sigmaVis(cohs==coh(n))];
        case 3
            % positive delta defined as ves to the left, vis to the right
            muVes = kves               * sind(hdg(n)-delta(n)/2) / 1000;
            muVis = kvis(cohs==coh(n)) * sind(hdg(n)+delta(n)/2) / 1000;

            % optimal weights (Drugo et al.) 
            wVes = sqrt( kves^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            wVis = sqrt( kvis(cohs==coh(n))^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            mu = wVes.*muVes + wVis.*muVis;

            % the DV is a sample from a dist with mean = weighted sum of
            % means. thus the variance is the weighted sum of variances
            % (error propagation formula):
            sigmaComb = sqrt(wVes.^2 .* sigmaVes^2 + wVis.^2 .* sigmaVis(cohs==coh(n))^2); % assume zero covariance
            s = [sigmaComb sigmaComb];
    end

    Mu = [mu,-mu]; % mean vector for 2D DV
     
    % convert correlation to covariance matrix
    V = diag(s)*S*diag(s);
    dv = [0 0; cumsum(mvnrnd(Mu,V,dur(n)-1))]; % bivariate normrnd
    % because Mu is signed according to heading (positive=right),
    % dv(:,1) corresponds to evidence favoring rightward, not evidence
    % favoring the correct decision (as in Kiani eqn. 3 and images_dtb)

    % decision outcome
    cRT1 = find(dv(:,1)>=B, 1);
    cRT2 = find(dv(:,2)>=B, 1);
    % the options are:
    % (1) only right accumulator hits bound,
    if ~isempty(cRT1) && isempty(cRT2)
        RT(n) = cRT1/1000 + Tnd;
        finalV(n) = dv(cRT1,2); % only 1 hit, so 2 is the loser
        hitBound(n) = 1;
        choice(n) = 1;
    % (2) only left accumulator hits bound, 
    elseif isempty(cRT1) && ~isempty(cRT2)
        RT(n) = cRT2/1000 + Tnd;
        finalV(n) = dv(cRT2,1); % only 2 hit, so 1 is the loser
        hitBound(n) = 1;
        choice(n) = -1;
    % (3) neither hits bound,
    elseif isempty(cRT1) && isempty(cRT2)
        RT(n) = dur(n)/1000 + Tnd;
            % this is interesting: which DV matters for confidence
            % if neither hits bound? for now take the (abs) maximum,
            % ie whichever was closest to hitting bound. (alternative
            % would be their average?)
        whichWon = dv(end,:)==max(dv(end,:));
        finalV(n) = dv(end,~whichWon); % the not-whichWon is the loser
        % % finalV(n) = mean(dvEnds);
        hitBound(n) = 0;
        a = [1 -1];
        choice(n) = a(whichWon);
    % (4) or both do
    else
        RT(n) = min([cRT1 cRT2])/1000 + Tnd;
        whichWon = [cRT1<=cRT2 cRT1>cRT2];
        finalV(n) = dv(min([cRT1 cRT2]),~whichWon); % the not-whichWon is the loser
        hitBound(n) = 1;
        a = [1 -1];
        choice(n) = a(whichWon);
    end
        
    % use map to look up log-odds that the motion is rightward
    diffV = abs((P.y+B)-finalV(n));
    diffT = abs(R.t-RT(n));
        
    thisV = find(diffV==min(diffV));
    thisT = find(diffT==min(diffT));
    logOddsCorr(n) = P.logOddsCorrMap(thisV(1), thisT(1));
    expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
    conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1
end

choice(choice==0) = sign(randn); % not needed under usual circumstances
choice(choice==1) = 2; choice(choice==-1) = 1; % 1=left, 2=right

% output var 'fit' gets the same values as data for the conditions, but the
% simulated trial outcomes for the observables!
fit.heading = data.heading;
fit.coherence = data.coherence;
fit.modality = data.modality;
fit.delta = data.delta;
fit.choice = choice;
fit.RT = RT;
fit.conf = conf;

fit.correct = (choice==2 & data.heading>0) | (choice==1 & data.heading<0);
data.correct = (data.choice==2 & data.heading>0) | (data.choice==1 & data.heading<0);


%% calculate error

% Do ibs_basic?

% % likelihood of rightward choice on each trial, under binomial assumptions
% Pr_model = 1 ./ (1 + exp(-2*k*B*sind(hdg)));

choiceD  = logical(data.choice-1);

Pr_model = nan(ntrials,1);
n = nan(length(mods),length(cohs),length(deltas),length(hdgs));

pRight = n;
RTmean_fit = n; RTmean_data = n; sigmaRT = n;
confMean_fit = n; confMean_data = n; sigmaConf = n;


for m = 1:length(mods)
for c = 1:length(cohs)
for d = 1:length(deltas)

    for h = 1:length(hdgs)
        J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d);
        
        n(m,c,d,h) = sum(J);
        nCor(m,c,d,h) = sum(J & data.correct); % use only correct trials for RT and conf fits
        
        pRight(m,c,d,h) = sum(J & fit.choice==2) / n(m,c,d,h); % 2 is rightward
        Pr_model(J) = pRight(m,c,d,h);
        
        RTmean_fit(m,c,d,h) = mean(fit.RT(J & data.correct));
        RTmean_data(m,c,d,h) = mean(data.RT(J & data.correct));
        sigmaRT(m,c,d,h) = std(data.RT(J & data.correct))/sqrt(nCor(m,c,d,h));
        
        confMean_fit(m,c,d,h) = mean(fit.conf(J & data.correct));
        confMean_data(m,c,d,h) = mean(data.conf(J & data.correct));
        sigmaConf(m,c,d,h) = std(data.conf(J & data.correct))/sqrt(nCor(m,c,d,h));
             
    end
end
end
end

% kluge to avoid log(0) issues
Pr_model(Pr_model==0) = min(Pr_model(Pr_model~=0)); 
Pr_model(Pr_model==1) = max(Pr_model(Pr_model~=1));
LL_choice = sum(log(Pr_model(choiceD))) + sum(log(1-Pr_model(~choiceD)));

% likelihood of mean RTs and mean confidence ratings for each condition, under Gaussian approximation

L_RT = 1./(sigmaRT*sqrt(2*pi)) .* exp(-(RTmean_fit - RTmean_data).^2 ./ (2*sigmaRT.^2));
L_RT(L_RT==0) = min(L_RT(L_RT~=0));
LL_RT = nansum(log(L_RT(:)));

L_conf = 1./(sigmaConf*sqrt(2*pi)) .* exp(-(confMean_fit - confMean_data).^2 ./ (2*sigmaConf.^2));
L_conf(L_conf==0) = min(L_conf(L_conf~=0));
LL_conf = nansum(log(L_conf(:)));

err = -(LL_choice + LL_RT + LL_conf);

%if call_num == 10, keyboard,end
%err = -LL_choice;

%% print progress report!
fprintf('\n\n\n****************************************\n');
fprintf('run %d\n', call_num);
fprintf('\tks= %g\n\tsigma= %g\n\tB= %g\n', param(1), param(2), param(3));
fprintf('err: %f\n', err);
if options.ploterr && strcmp(options.fitMethod,'fms')
    figure(options.fh); hold on
    plot(call_num, err, '.','MarkerSize',12,'color','k');
    drawnow;
end
call_num = call_num + 1;



end


