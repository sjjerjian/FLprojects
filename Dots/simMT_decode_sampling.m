% very simple sampling-based decoding

dirAxis = prefDirs; % here it's only defined at the neurons' pref dir, not for all possible dirs

clr = {'k-','b-','g-','c-','m-','r-'};
Posterior = nan(nTrials,length(dirAxis));
choice = nan(nTrials,1);
pdw = nan(nTrials,1);
vecSum = nan(nTrials,1);

% for fitting of circular Gaussian (Von Mises) to likelihood/posterior
options = optimset('MaxFunEvals', 1e9, 'Algorithm','interior-point','LargeScale','off','Display','off');
% options = optimset('PlotFcns',@optimplotfval);
vonmises_fcn = @(b,x) b(1) * exp(b(2)*cosd(x-b(3))) / (2*pi*besseli(0,b(2)));
err_fcn = @(param,X,data) sum((vonmises_fcn(param,X)-data).^2); % sum sq. error
kappaGuess = 2; % inverse variance term
betaAll = nan(nTrials,3); % matrix for storing the fitted params
errAll = nan(nTrials,1);

gamma = 0.4; % threshold on inverse variance for betting high (arbitrary)

tic
for t = 1:nTrials
            
    % simplest version of sampling I can think of is that every spike is a
    % sample from the posterior at the value defined by the cell's pref dir
%     P = zeros(1,length(prefDirs));
%     for s = 1:size(R,3)
%         P = P + R(t,:,s);
%         % ...
%     end
    % this can give you P at any given time point. but in the initial case
    % we'll just wait for stimulus offset and use the summed spikes as the
    % posterior. Thus we only need to sum the spikes, and we simply use the
    % familiar 'hill' of activity, fitted, to get the posterior width (?!)
    % This can't be right...
    hill = nansum(squeeze(R(t,:,:))'); %#ok<UDIM>
    P = hill/sum(hill);
    
    % I guess the alternative is that each spike counts as a sample of a
    % multivariate posterior which is updated by sampling from its
    % conditional distributions, as in Gibbs sampling. Each cell defines a
    % random variable and the full posterior is the joint posterior across
    % all cells/variables. This sounds more correct but not sure how to do
    % it.
    
    % fit von mises to get width for confidence, and plot a few examples
    beta_guess = [max(P) kappaGuess dir(t)];
        % % test guess
        % figure; plot(dirAxis,P,'b-',dirAxis,vonmises_fcn(beta_guess,dirAxis),'r--')
    
    [beta,err] = fminsearch(@(j) err_fcn(j,dirAxis,P), beta_guess, options);    
        
%         % temp: check quality of fit
%         % (tends to be poor with coh near zero -- why???)
%         if t<11
%             figure(t); plot(dirAxis,P,'b-',dirAxis,vonmises_fcn(beta,dirAxis),'r--')
%         end
        
    % choice based on MAP estimate (min sq error cost fn)
    if sign(cosd(beta(3)))==0
        choice(t) = sign(randn);
    else
        choice(t) = sign(cosd(beta(3)));
    end

%     % alternate metric for confidence: magnitude of vector sum (not good)
%     [x,y] = pol2cart(dirAxis*pi/180,P);
%     vecSum(t) = sqrt(sum(x)^2 + sum(y)^2);
    
    betaAll(t,:) = beta;
    errAll(t) = err;
    Posterior(t,:) = P;

end
toc

P = Posterior;

choice = (choice+1)/2; % convert to 0/1

RT = nan(size(choice)); % no mechanism for RT in PPC

% confidence based on width of posterior;
% convert to wager by comparison to arbitrary threshold
pdw(betaAll(:,2)>=gamma) = 1;
pdw(betaAll(:,2)<gamma) = 0;



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
    % why doesn't it vary with dur?



