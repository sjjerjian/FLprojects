% fiddle_goris

clear; close all
load temp.mat

% histogram of spike counts conditioned on a stimulus
coh = 0;
m=1; c=1;
nbins = 20;

I = dataCell{m}.Exp.openEvents.coherence==coh;
count = spCount_exp{m}(c,I);
rate = spRate_exp{m}(c,I);
figure(1); hist(rate,nbins); % can't use count, because a given distribution is defined by a given dt, and ours varies across trials
[n,~] = hist(count,nbins);

%% poisson dist
dt = RT(I)'/1000; 
mu = mean(count);
    % simplify for now: use rate instead of count and and call dt=1
muRate = mean(rate);
N = 1:max(count)+3;
pN = (muRate.^N)./factorial(N) * exp(-muRate);
hold on; plot(N,pN*max(n)/max(pN),'b');

% now fit it
poiss_err = @(mu,N,dt) -sum(log(((mu*dt).^N)./factorial(N) .* exp(-mu*dt))); % error is negative log likelihood of observing the data given the param
% poiss_err(muRate,count,dt) % sanity check: should be near minimum with mu=muRate, and go up in both directions (mu>muRate and mu<muRate)
muFit = fminsearch(@(mu) poiss_err(mu,count,dt), muRate);
pNfit = (muFit.^N)./factorial(N) * exp(-muFit);
    % ^still not sure how to compute the smooth curve, because dt can't be 
    % same length as N (so I just omit it, which makes no sense)...
    % but oh well, moving on
plot(N,pNfit*max(n)/max(pNfit),'r');

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
sigmaG = 0.2;
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
x = [0.2 muRate];
X = fminsearch(@(x) modPoiss_err(x,rate), x); % hooray!

pNmod = modPoiss(X(1),X(2),N);
hold on; plot(N,pNmod*max(n)/max(pNmod),'m');


%% okay now we can check whether sigmaG varies with, say, coherence

m=1; c=1;
nbins = 20;
dir = 0;

ucoh = unique(dataCell{m}.Exp.openEvents.coherence);
for C = 1:length(ucoh)

    I = dataCell{m}.Exp.openEvents.coherence==ucoh(C) & dataCell{m}.Exp.openEvents.direction==dir; % check other dir
    rate = spRate_exp{m}(c,I);
    muRate = mean(rate);
    x = [0.2 muRate];
    X_all(C,:) = fminsearch(@(x) modPoiss_err(x,rate), x);

    % now check PDW
    I = dataCell{m}.Exp.openEvents.coherence==ucoh(C) & dataCell{m}.Exp.openEvents.direction==dir & dataCell{m}.Exp.openEvents.pdw==1;
    rate = spRate_exp{m}(c,I);
    muRate = mean(rate);
    x = [0.2 muRate];
    X_high(C,:) = fminsearch(@(x) modPoiss_err(x,rate), x);

    I = dataCell{m}.Exp.openEvents.coherence==ucoh(C) & dataCell{m}.Exp.openEvents.direction==dir & dataCell{m}.Exp.openEvents.pdw==0;
    rate = spRate_exp{m}(c,I);
    muRate = mean(rate);
    x = [0.2 muRate];
    X_low(C,:) = fminsearch(@(x) modPoiss_err(x,rate), x);
    
end

figure;plot(ucoh,X_all(:,1),'b-o',ucoh,X_low(:,1),'r-o',ucoh,X_high(:,1),'g-o');
legend('all','bet-low','bet-high')
















