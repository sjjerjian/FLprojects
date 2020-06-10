% very simple PPC/likelihood-based decoding

dirAxis = 0:360;
clr = {'k-','b-','g-','c-','m-','r-'};
P = nan(nTrials,length(dirAxis));
choice = nan(nTrials,1);
pdw = nan(nTrials,1);
vecSum = nan(nTrials,1);

% for fitting of circular Gaussian (Von Mises) to likelihood/posterior
options = optimset('MaxFunEvals', 1e9, 'Algorithm','interior-point','LargeScale','off','Display','off');
% options = optimset('PlotFcns',@optimplotfval);
vonmises_fcn = @(b,x) b(1) * exp(b(2)*cosd(x-b(3))) / (2*pi*besseli(0,b(2)));
err_fcn = @(param,X,data) sum((vonmises_fcn(param,X)-data).^2); % sum sq. error
kappaGuess = 5; % inverse variance term
betaAll = nan(nTrials,3); % matrix for storing the fitted params
errAll = nan(nTrials,1);

gamma = 10; % threshold on inverse variance for betting high (arbitrary)

% The model assumes a readout mechanism that 'knows' the tuning curves.
% The idea is to compute the likelihood function (P(r|S) as a function
% of S) by, in effect, comparing the response on a given trial with the
% expected responses given by the tuning curves (see likelihood_tutorial.m).
% Since tuning{} is specified in spikes/s, what we need is a vector of
% spike rates (not counts) for each trial:
Rates = nan(nTrials,nNeurons);
for t = 1:nTrials
    win = latency:dur(t)+latency;
    Rates(t,:) = nansum(R(t,:,win),3) / (dur(t)/1000);
end

h = cell(length(poscohs),1);
for c = 1:length(poscohs)
    % 'kernel' (h) for independent poisson is log of the tuning curves
    h{c} = log(tuning{c}/10);    
                    %^ divide by 10 to keep the exp() below from blowing up
                    % to inf; will compensate by dividing r by 10 also
                       
        % or, take into account covariance matrix [work in progress; see Ma
        % et al. 2006 Eq. 5 (I still don't really understand it)
%     I = abs(coh)==poscohs(c);
%     h2{c} = inv(cov(Rates(I,:))) * diff(tuning{c},1,2);
end


tic
for t = 1:nTrials
        
    c = poscohs == abs(coh(t));
    r = Rates(t,:);
    L = exp(transpose(h{c}) * (r'/10)); % (divide by 10 here too; hope this is okay!)
    
    % normalized likelihood; this is the posterior, assuming a flat prior
    Lnorm = L/sum(L);

    % fit von mises to get width for confidence, and plot a few examples
    beta_guess = [max(Lnorm) kappaGuess dir(t)]; 
        % % test guess
        % figure; plot(dirAxis,Lnorm,'b-',dirAxis,vonmises_fcn(beta_guess,dirAxis),'r--')

    [beta,err] = fminsearch(@(j) err_fcn(j,dirAxis,Lnorm'), beta_guess, options);    
        
%         % temp: check quality of fit
%         % (tends to be poor with coh near zero -- why???)
%         if t<11
%             figure(t); plot(dirAxis,Lnorm,'b-',dirAxis,vonmises_fcn(beta,dirAxis),'r--')
%         end
        
    % choice based on log likelihood ratio for the two alternatives
    choice(t) = sign(log(L(dirAxis==0)/L(dirAxis==180)));
    %[would MAP estimate (+cost function) be different? I don't think so..]

%     % alternate metric for confidence: magnitude of vector sum (not good)
%     [x,y] = pol2cart(dirAxis*pi/180,Lnorm');
%     vecSum(t) = sqrt(sum(x)^2 + sum(y)^2);
    
    betaAll(t,:) = beta;
    errAll(t) = err;
    P(t,:) = Lnorm;

end
toc

% % trying to diagnose/improve quality of fits
% nanmean(errAll)
% nanmedian(errAll)

choice = (choice+1)/2; % convert to 0/1

RT = nan(size(choice)); % no mechanism for RT in PPC

% confidence based on width of normalized likelihood (posterior);
% convert to wager by comparison to arbitrary threshold
pdw(betaAll(:,2)>=gamma) = 1;
pdw(betaAll(:,2)<gamma) = 0;

% % alternative (not a good one)
% pdw(vecSum>=0.92) = 1;
% pdw(vecSum<0.92) = 0;



%% some sanity checks

% posterior as a function of coherence
for c = 1:length(poscohs)
    I = coh==-poscohs(c); % use only 180 dir (neg coh) for easier visualization
    figure(14); plot(mean(P(I,:)),clr{c}); hold on;
end
title('average posterior for each unsigned coh (true dir=180)')
xlabel('dir'); ylabel('P(dir|R)'); xlim([0 360]);

% posterior as a function of duration
D = [0 quantile(dur,5) inf];
for d = 1:length(D)-1
    I = dur>=D(d) & dur<D(d+1) & dir==180;
    figure(15); plot(mean(P(I,:)),clr{d}); hold on;
end
title('average posterior for quantiles of duration (true dir=180)')
xlabel('dir'); ylabel('P(dir|R)'); xlim([0 360]);



