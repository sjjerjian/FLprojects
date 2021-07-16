% fiddle_goris

clear; close all
load temp2.mat

T = 1; U = 0;
for n = 1:length(dataCell)
    n
    for c = 1:nCh(n)
%         if max(spCount_exp(n,c,:))<147 % 170? [max to avoid infs, shouldn't exclude many cells];
        if max(spRate_exp(n,c,:))<147 % 170? [max to avoid infs, shouldn't exclude many cells];
        U = U+1;
            for t = 1:length(dataCell{n}.Exp.openEvents.coherence)
                if ~isnan(spRate_exp(n,c,t)) && ~isnan(spCount_exp(n,c,t))
                    sprate(T,1) = spRate_exp(n,c,t); 
                    spcount(T,1) = spCount_exp(n,c,t);
        % %             sprate_check(T,1) = spcount(T,1)/(RT(T,1)/1000); 

                    unit(T,1) = U; 
                    coh(T,1) = dataCell{n}.Exp.openEvents.coherence(t);
                    dir(T,1) = dataCell{n}.Exp.openEvents.direction(t); 
                    choice(T,1) = dataCell{n}.Exp.openEvents.choice(t); 
                    pdw(T,1) = dataCell{n}.Exp.openEvents.pdw(t);            
                    mStart = dataCell{n}.Exp.openEvents.motionStart(t);            
                    mEnd = dataCell{n}.Exp.openEvents.motionEnd(t);
                    RT(T,1) = round((mEnd-mStart)*1000);        
                    T = T+1;
                end
            end   
        end
    end
end

pref = nan(size(unit));
normrate = nan(size(sprate));
for U = 1:max(unit)
    U
    I = unit==U & dir==0 & coh>0.25;
    J = unit==U & dir==180 & coh>0.25;
    if nanmean(sprate(I)) > nanmean(sprate(J))
        pref(unit==U) = 0;
    else
        pref(unit==U) = 180;
    end

    normrate(unit==U) = sprate(unit==U)/max(sprate(unit==U));
    
end

%% poisson dist

% will be plotting a histogram of spike counts (with nbins)
% conditioned on a particular stimulus strength (whichcoh)
whichcoh = 0;
nbins = 20;

I = coh==whichcoh; % dir doesn't matter if whichcoh==0 
count = spcount(I);
rate = sprate(I);
figure(1); hist(rate,nbins); % can't use count, because a given distribution is defined by a given dt, and ours varies across trials
[n,~] = hist(count,nbins);

dt = RT(I)'/1000; 
mu = mean(count);
    % simplify for now: use rate instead of count and and call dt=1
muRate = mean(rate);
N = 1:max(count);
pN = (muRate.^N)./factorial(N) * exp(-muRate);
hold on; plot(N,pN*max(n)/max(pN),'b');


% % now fit it . -**not working for some reason
% poiss_err = @(mu,N,dt) -sum(log(((mu*dt).^N)./factorial(N) .* exp(-mu*dt))); % error is negative log likelihood of observing the data given the param
% % poiss_err(muRate,count,dt) % sanity check: should be near minimum with mu=muRate, and go up in both directions (mu>muRate and mu<muRate)
% muFit = fminsearch(@(mu) poiss_err(mu,count,dt), muRate);
% pNfit = (muFit.^N)./factorial(N) * exp(-muFit);
%     % ^still not sure how to compute the smooth curve, because dt can't be 
%     % same length as N (so I just omit it, which makes no sense)...
%     % but oh well, moving on
% plot(N,pNfit*max(n)/max(pNfit),'r');


% alt:
% poiss_err = @(muR,N) -sum(log((muR.^N)./factorial(N) .* exp(-muR)));
% muFit = fminsearch(@(muR) poiss_err(muR,count), muRate);    
% pNfit = (muFit.^N)./factorial(N) * exp(-muFit);
% hold on; plot(N,pNfit*max(n)/max(pNfit),'m');


%% modulated poisson (Goris et al. 2014)

modPoiss = @(sigmaG,mu,N) gamma(N+(1/sigmaG^2)) ./ (gamma(N+1).*gamma(1/sigmaG^2)) .* ...
                             (sigmaG^2 * mu).^N ./ ((sigmaG^2 * mu + 1).^(N+1/sigmaG^2));
        % OR, w delta-t:
% modPoiss = @(sigmaG,mu,dt,N) gamma(N+(1/sigmaG^2)) ./ (gamma(N+1).*gamma(1/sigmaG^2)) .* ...
%                              (sigmaG^2*mu*dt).^N ./ ((sigmaG^2*mu*dt + 1).^(N+1/sigmaG^2));

% try hand-tuning to find a rough guess for sigmaG
sigmaG = 0.9;

pNmod = modPoiss(sigmaG,muRate,N);
% OR
% pNmod = modPoiss(sigmaG,mu,mean(dt),N);

figure(2); hist(rate,nbins);
hold on; plot(N,pNmod*max(n)/max(pNmod),'g');


% now fit it

% modPoiss_err = @(sigmaG,mu,N) -sum(log(modPoiss(sigmaG,mu,N)));
% % modPoiss_err(0.2,muRate,rate) % sanity check: should be near minimum with mu=muRate and sigmaG~=0.2, and go up as you deviate from those
% 
% % param = fminsearch(@(sigmaG,muRate) modPoiss_err(sigmaG,muRate,rate), [0.2 muRate]);
%     % these two below work fine, so why the error msg with fminsearch?
%     modPoiss_err(sigmaG,muRate,rate)
%     modPoiss(sigmaG,muRate,rate)
% 
% sigGfit = fminsearch(@(sigmaG) modPoiss_err(sigmaG,muRate,rate), 0.2); % this works... 

% so let's try a vector of params instead of passing them in individually:
modPoiss_err = @(param,N) -sum(log(modPoiss(param(1),param(2),N)));
guess = [0.9 muRate];

% % debugging
% modPoiss_err(guess,rate)
% log(modPoiss(guess(1),guess(2),N))

X = fminsearch(@(x) modPoiss_err(x,rate), guess); % hooray!

pNmod = modPoiss(X(1),X(2),N);
hold on; plot(N,pNmod*max(n)/max(pNmod),'m');


%% okay now we can check whether sigmaG varies with coherence and PDW!

ucoh = unique(coh);
for C = 1:length(ucoh)

    I = coh==ucoh(C) & dir==pref & ~isnan(sprate);
    rate = sprate(I);
    muRate = mean(rate);
    guess = [0.9 muRate];
    X_all(C,:) = fminsearch(@(x) modPoiss_err(x,rate), guess);

    % now check PDW
    I = coh==ucoh(C) & dir==pref & pdw==1 & ~isnan(sprate);
    rate = sprate(I);
    muRate = mean(rate);
    guess = [0.9 muRate];
    X_high(C,:) = fminsearch(@(x) modPoiss_err(x,rate), guess);

    I = coh==ucoh(C) & dir==pref & pdw==0 & ~isnan(sprate);
    rate = sprate(I);
    muRate = mean(rate);
    guess = [0.9 muRate];
    X_low(C,:) = fminsearch(@(x) modPoiss_err(x,rate), guess);
    
end

figure;plot(ucoh,X_all(:,1),'b-o',ucoh,X_low(:,1),'r-o',ucoh,X_high(:,1),'g-o');
legend('all','bet-low','bet-high'); title('pref dir trials');


for C = 1:length(ucoh)

    I = coh==ucoh(C) & dir~=pref & ~isnan(sprate);
    rate = sprate(I);
    muRate = mean(rate);
    guess = [0.9 muRate];
    X_all(C,:) = fminsearch(@(x) modPoiss_err(x,rate), guess);

    % now check PDW
    I = coh==ucoh(C) & dir~=pref & pdw==1 & ~isnan(sprate);
    rate = sprate(I);
    muRate = mean(rate);
    guess = [0.9 muRate];
    X_high(C,:) = fminsearch(@(x) modPoiss_err(x,rate), guess);

    I = coh==ucoh(C) & dir~=pref & pdw==0 & ~isnan(sprate);
    rate = sprate(I);
    muRate = mean(rate);
    guess = [0.9 muRate];
    X_low(C,:) = fminsearch(@(x) modPoiss_err(x,rate), guess);
    
end

figure;plot(ucoh,X_all(:,1),'b-o',ucoh,X_low(:,1),'r-o',ucoh,X_high(:,1),'g-o');
legend('all','bet-low','bet-high'); title('null dir trials');



%% but wait, resps are not normalized but pooled.
% either need to normalize, or compute cell by cell and average.


% use normalized rate
sprate = normrate*50;

ucoh = unique(coh);
for C = 1:length(ucoh)

    I = coh==ucoh(C) & dir==pref & ~isnan(sprate);
    rate = sprate(I);
    muRate = mean(rate);
    guess = [0.9 muRate];
    X_all(C,:) = fminsearch(@(x) modPoiss_err(x,rate), guess);

    % now check PDW
    I = coh==ucoh(C) & dir==pref & pdw==1 & ~isnan(sprate);
    rate = sprate(I);
    muRate = mean(rate);
    guess = [0.9 muRate];
    X_high(C,:) = fminsearch(@(x) modPoiss_err(x,rate), guess);

    I = coh==ucoh(C) & dir==pref & pdw==0 & ~isnan(sprate);
    rate = sprate(I);
    muRate = mean(rate);
    guess = [0.9 muRate];
    X_low(C,:) = fminsearch(@(x) modPoiss_err(x,rate), guess);
    
end

figure;plot(ucoh,X_all(:,1),'b-o',ucoh,X_low(:,1),'r-o',ucoh,X_high(:,1),'g-o');
legend('all','bet-low','bet-high'); title('pref dir trials');


for C = 1:length(ucoh)

    I = coh==ucoh(C) & dir~=pref & ~isnan(sprate);
    rate = sprate(I);
    muRate = mean(rate);
    guess = [0.9 muRate];
    X_all(C,:) = fminsearch(@(x) modPoiss_err(x,rate), guess);

    % now check PDW
    I = coh==ucoh(C) & dir~=pref & pdw==1 & ~isnan(sprate);
    rate = sprate(I);
    muRate = mean(rate);
    guess = [0.9 muRate];
    X_high(C,:) = fminsearch(@(x) modPoiss_err(x,rate), guess);

    I = coh==ucoh(C) & dir~=pref & pdw==0 & ~isnan(sprate);
    rate = sprate(I);
    muRate = mean(rate);
    guess = [0.9 muRate];
    X_low(C,:) = fminsearch(@(x) modPoiss_err(x,rate), guess);
    
end

figure;plot(ucoh,X_all(:,1),'b-o',ucoh,X_low(:,1),'r-o',ucoh,X_high(:,1),'g-o');
legend('all','bet-low','bet-high'); title('null dir trials');











