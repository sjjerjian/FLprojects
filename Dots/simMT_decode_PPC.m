% very simple PPC/likelihood-based decoding

dirAxis = 0:360;
clr = {'k-','b-','g-','c-','m-','r-'};
P = nan(nTrials,length(dirAxis));
P_Kernel = nan(nTrials,length(dirAxis)-1);
choice = nan(nTrials,1);
choice_Kernel = nan(nTrials,1); %For use with Kernel that takes into account correlation
choice_Kernel_v2 = nan(nTrials,1);
pdw = nan(nTrials,1);
var_Kernel = nan(nTrials,1); %For use with Kernel that takes into account correlation
vecSum = nan(nTrials,1);
%space for Time based PPC
choice_Kernel_Time = nan(nTrials,1);
time_Kernel_Time = nan(nTrials,1);
var_Kernel_Time = nan(nTrials,1);
%mvl: A lot of new variables simply for storaging purposes


%mvl: von mises original paramters
dirAxis = 0:359;
K = 3; % inverse variance term (will scale with coh)... mvl: Read Primer paper (might be in there, how to get K)
ampl = 60;  % actual peak FR will be about half of this, at 51.2% coh
offset = 0;
prefDirs = linspace(dirAxis(1),dirAxis(end),nNeurons); %mvl: Pref direction for all the neurons (One neuron for each direction...)
for c = 1:length(poscohs) %mvl: Tuning curve for Each Neuron... And For Each Coherence.
    vonmises_Der{c} = zeros(nNeurons,length(dirAxis)); %mvl: Number of Neurons x Each Direction (Stimulus)... Each Neurons tuning curve. 
    for n = 1:nNeurons
        % von Mises tuning functions (captures asymmetry in pref/left response to coh)
        k = K * poscohs(c); % scale inverse variance by coh
        vonmises_Der{c}(n,:) = ampl .* -k*sind(dirAxis - prefDirs(n)) .* exp(k*cosd(dirAxis-prefDirs(n))) / (2*pi*besseli(0,k)) + offset; %von mises derivative, for kernel function (derived by hand)
                 %^ rows are neurons, columns are dirs   '
    end
    if 0
        plot(vonmises_Der{c}(round(nNeurons/2),:)); axis tight; hold on;  %look at the derivative function
    end
end
if 0
    figure;
    plot(dirAxis, vonmises_Der{6}(1:round(nNeurons/10):nNeurons, :)); axis tight;
    xlabel('Direction')
end
%% mvl: found the derivatives of the runing curve, now find the kernel
% first reconstruct covariance or reuse the correlation matrix (is reconstructing then
% you need to use the given data, but for now just use the given
% correlation). CAVEAT: This method of decoding assumes the neurons have an
% understanding of their own correlation. No?

% first create corr mat
maxCorr = 0.2;
    % pairwise correlation varies smoothly from maxCorr to zero based on
    % tuning similarity (pref dir proximity)
fano = 1.8; % fano factor = variance/mean ratio (forced to be 1 for poisson,
            % but for real neurons it's often 1.5-1.8)
c=linspace(maxCorr,0.0,nNeurons/2); 
c=[c fliplr(c)];
Cor=toeplitz(c); % aka diagonal-constant matrix, shortcut for making the kind of corr mat we want
Cor(logical(eye(size(Cor)))) = 1; % set main diagonal to be 1
tmpCor = Cor; 
tmpCor(tmpCor>.9)=NaN; % for better color range when plotting

%We will try using true correlation first, if this doesnt give good answers then
%reconstruct the covariance matrix (In the Beck 2008 paper, a continuation
%of the Ma 2006 one, it is clear that each Coherence warrents a different
%Covariance and Tuning curve, since both are dependent on the lvl of coh.
%However that equation can be modified h'(s) = 1/(g(c)*cov(s)) *
%(g(c)*f'(s)), so that the underlying correlation or covariance is
%independent of a gain (coherence in this case)... I think this is right.
%But again in my case just using True Correlation (the one chris
%constructed originally)
dirAxis = 0:359;
for c = 1:length(poscohs)
    kernelH_Der{c} = Cor\vonmises_Der{c}; %mvl: h'(s) = inv(Sigma)*f'(s)
    kernelH{c} = cumsum(kernelH_Der{c}, 2); %mvl: intergral-> integrate through Stimulus/row (there will errors since this is numerical estimate), 
    %We have to reduce the kernel because it is increasing numbers too high
    %and giving nan and inf, once again this might be mention in Wei Ji Ma
    %2006 paper, as removing the Gain, or making the Kernel's estimation by
    %Covariance and Fitting Function independent of Gain
    %Another thing, ideally all kernels should have the same max amplitude and min.
    %However, this is not the case so we will manually fix this by having
    %the max amplitude be the same for all kernels (this could be wrong)
    kernelH{c} = (kernelH{c} + (10 - max(kernelH{c},[],2)))./500; %mvl: Instead, offset the Kernel so max amp matches and then scale to reduce load on computer (normally there is some weird mislaignment between Neurons, some have bigger max amps then other, why?) I think it has to do with the covariance matrix structure and its inverse.
    if 1
        figure(100);
        plot(dirAxis, kernelH{c}(round(nNeurons/2), :)); hold on; 
        figure(101);
        plot(dirAxis, log(tuning{c}(round(nNeurons/2),:))); hold on;
    end
end
hold off
% plot kernel to see if it looks anything like the Log(Vonmises) that
% independent poisson gives (shapes look similar but values are vastly
% different, because i divided by 500. When you add it it looks better).
if 1
    figure
    plot(dirAxis, log(tuning{6}(round(nNeurons/2),:))); hold on;
    plot(dirAxis, tuning{6}(round(nNeurons/2),:)); hold off;
    title('Kernel with NO Correlation, with Log')
    xlabel('Direction')

    figure
    plot(dirAxis, kernelH{1}(round(nNeurons/2), :).*500)
    title('Kernel with Correlation')
    xlabel('Direction')
    
end

%mvl: Another reason I think I limit and made every kernel to have the same
%maximum amplitude was because we want each's neurons Kernel to have
%similar effects. If one Neuron's kernel is bigger in magnitude (amplitude
%of kernel) then that neuron will always contribute more and we dont want
%that).
%%
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
    Rates(t,:) = nansum(R(t,:,win),3) / (dur(t)/1000); %mvl: convert spikes into Rates of second (spikes/second)
end

h = cell(length(poscohs),1);
for c = 1:length(poscohs)
    % 'kernel' (h) for independent poisson is log of the tuning curves
    h{c} = log(tuning{c}/10);  %mvl: Log of the Mean FR for each Neuron for Each direction (by Coherence), I believe this is used to weight things  
                    %^ divide by 10 to keep the exp() below from blowing up
                    % to inf; will compensate by dividing r by 10 also
                       
        % or, take into account covariance matrix [work in progress; see Ma
        % et al. 2006 Eq. 5 (I still don't really understand it)
%     I = abs(coh)==poscohs(c);
%     h2{c} = inv(cov(Rates(I,:))) * diff(tuning{c},1,2);
end


tic
%mvl: Beck 2008
% We first want to bin the Responses R by 50ms
binT = 50;
binR = []; %bin Rates
binS = []; %bin spikes
for i = 1:floor(size(R,3)/50)
    tempBinR = nansum(R(:,:, (binT*(i-1))+1: binT*i), 3);
    binS(:,:,i) = tempBinR;
    binR(:,:,i) = tempBinR./(50/1000); %use rates, see how much of a difference that makes
end
%Create a matrix or cell that can contain all the posteriors for LIP, that
%takes into account the summation for each time bin
suffStat_Time = cell(1,nTrials);
Norm_suffStat_Time = cell(1,nTrials);
bound_Kernel_Time = .1;

dirAxis = 0:360;
for t = 1:nTrials
        
    c = poscohs == abs(coh(t));
    r = Rates(t,:); %mvl: neuron response (for this trial) 
    L = exp(transpose(h{c}) * (r'/10)); % (divide by 10 here too; hope this is okay!), mvl: weighting occuring here
    suffStat = exp(kernelH{c}' * (r'/10)); % mvl: exp(h'(s)*r), it might be blowing up to inf because there is the need to remove the Gain (In supplementary material, on how to go about this),
    %mvl: try using spike counts instead of spike Rates see how that changes this
    %crude method for now.
    %mvl: if I remembe correctly scaling the KernelH helped.
    
    %mvl: Beck 2008
    hitbound = 0;
    for ii = 1:size(binR,3)
        r = sum(binS(t,:,1:ii), 3); % Sum(r_MT(t_n))]
        %Mvl: dont complete understand how to work the null, keep this for
        %now since it works
        findNull = kernelH_Der{c}(:,dir(t)+1)'; %Find the null for this vector (for "global inhibition", I dont think i got this to work properly)
        N_KernelH_Der = null(findNull); %mvl: Also construct the nullspace h'(s)*n = 0, null(A) gives a list of null vectors... pick one???
        index_0 = find(findNull*N_KernelH_Der==0, 1); %mvl: Find the null that gives exactly 0
        suffStat_Time{t}(ii,:) = exp(kernelH{c}' * (r' - N_KernelH_Der(:,index_0))); %mvl: subtract by the null, its suppose to reduce saturation 
        Norm_suffStat_Time{t}(ii,:) = suffStat_Time{t}(ii,:)./sum(suffStat_Time{t}(ii,:),2); %mvl: normalize the distribution 
        if hitbound == 0
            if any(max(Norm_suffStat_Time{t}(ii,:)) >= bound_Kernel_Time)
                [~, tempDeg] = max(Norm_suffStat_Time{t}(ii,:)); %Find angle that his bound, that is angle in the distribution hits the bound mark the time, use max to avoid two hits or more. Pick the biggest one
                choice_Kernel_Time(t) = tempDeg - 1; %Best Stimulus approximation
                time_Kernel_Time(t) = ii;
                var_Kernel_Time(t) = 1/max(Norm_suffStat_Time{t}(ii,:)); %certainty for now is the inverse of the max (from 2006 Ma paper)
                hitbound = 1; %never come back to these if statements)
            elseif ii == size(binR,3) %we have come to the end and nothing has been hit then just choose the direction that contains the max value
                [~, choice_Kernel_Time(t)] = max(Norm_suffStat_Time{t}(ii,:)); %Find angle that his bound, that is angle in the distribution hits the bound mark the time, use max to avoid two hits or more. Pick the biggest one
                time_Kernel_Time(t) = ii;
                var_Kernel_Time(t) = 1/max(Norm_suffStat_Time{t}(ii,:));
            end           
        end
    end
    
    % normalized likelihood; this is the posterior, assuming a flat prior
    Lnorm = L/sum(L);
    ssNorm = suffStat/sum(suffStat); %mvl: possibly the eta in the equation. What happens to the phi(r,g) though? All i know is that it can (i guess) be ignored since it is not dependent on the stimulus
    

    % fit von mises to get width for confidence, and plot a few examples
    beta_guess = [max(Lnorm) kappaGuess dir(t)]; 
        % % test guess
        % figure; plot(dirAxis,Lnorm,'b-',dirAxis,vonmises_fcn(beta_guess,dirAxis),'r--')

    [beta,err] = fminsearch(@(j) err_fcn(j,dirAxis,Lnorm'), beta_guess, options);
    
    %mvl: instead of fitting von mises just fit a gaussian distribution
    %(will give mean and variance). This is so because in the paper I
    %believe it is stated that a combination of poisson distibutions (in
    %the factor format) gives rise to a gaussian. I believe it works for
    %mean of direction, but working with a circle makes it hard to computer
    %the second moment.
    firstMoment_SS = ssNorm' * [cosd(0:359)' sind(0:359)'];%dirAXis
    meanSS = atan2(firstMoment_SS(2), firstMoment_SS(1)); %simply equation for finding angle given Sin and Cos
    meanSS = rad2deg(mod(meanSS,2*pi)); %trick to convert to radiant in all 4 quadrants, then to angle
    secondMoment_SS = ssNorm' * [cos(((0:359).*(pi./180)).^2)' sin(((0:359).*(pi./180)).^2)'];
    secondMoment_SS = atan2(secondMoment_SS(2), secondMoment_SS(1));
    secondMoment_SS = rad2deg(mod(secondMoment_SS,2*pi));
    varSS = secondMoment_SS - meanSS.^2;
    % mod(atan2(-.9397,.3420),2*pi);
    
    %We can also make a choice using sign of Log like Chris did in the
    %bottom
    choice_Kernel_v2(t) =  sign(log(ssNorm(dirAxis==0)/ssNorm(dirAxis==180)));
    % sum(choice_Kernel_v2' == cosd(dir)) %Checking results purposes
    
    
%         % temp: check quality of fit
%         % (tends to be poor with coh near zero -- why???)
%         if t<11
%             figure(t); plot(dirAxis,Lnorm,'b-',dirAxis,vonmises_fcn(beta,dirAxis),'r--')
%         end
        
    % choice based on log likelihood ratio for the two alternatives
    choice(t) = sign(log(L(dirAxis==0)/L(dirAxis==180)));
    choice_Kernel(t) = meanSS;
    %sum(sign(cosd(choice_Kernel)) == cosd(dir)) %checking purposes
    var_Kernel(t) = 1/max(ssNorm); %mvl: Variance of 'posterior' distribution
    %[would MAP estimate (+cost function) be different? I don't think so..]

%     % alternate metric for confidence: magnitude of vector sum (not good)
%     [x,y] = pol2cart(dirAxis*pi/180,Lnorm');
%     vecSum(t) = sqrt(sum(x)^2 + sum(y)^2);
    
    betaAll(t,:) = beta;
    errAll(t) = err;
    P(t,:) = Lnorm;
    
    %mvl:Safe 'posterior' for every trial
    P_Kernel(t,:) = ssNorm;
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


%mvl: certainty measurement based of Norm Variance Measurement of
%distribution
% norm_var_kernel = normalize(var_Kernel);
% certainty = nan(size(norm_var_kernel));
% gamma_Kernel = mean([max(norm_var_kernel) min(norm_var_kernel)]);
% certainty(norm_var_kernel < gamma_Kernel) = 1;
% certainty(norm_var_kernel > gamma_Kernel) = 0;

%mvl: The above doesnt work that well, need to understand how to get second
%moment from circular distribution. For now use certainty = 1/Gain. Where,
%for now, Gain = Max of ssNorm or Posterior (crude amplitude calculation).
certainty = nan(size(var_Kernel));
gamma_Kernel = mean([max(var_Kernel) min(var_Kernel)]);
certainty(var_Kernel < gamma_Kernel) = 1;
certainty(var_Kernel > gamma_Kernel) = 0;
% Graph both Choice and Variance based on MVL Computations
avgChoice_Kernel = nan(size(cohs));
avgCertainty_Kernel = nan(size(cohs));
avgChoice_Kernel_Time = nan(size(cohs));
avgCertainty_Kernel_Time = nan(size(cohs));
avgRT_Kernel_Time = nan(size(cohs));
        
if 1
    for i = 1:length(cohs)
        avgChoice_Kernel(i) = mean((sign(cosd(choice_Kernel(coh == cohs(i)))) + 1)./2);
        avgCertainty_Kernel(i) = mean(certainty(coh == cohs(i)));
        
        avgChoice_Kernel_Time(i) = mean((sign(cosd(choice_Kernel_Time(coh == cohs(i)))) + 1)./2);
        avgCertainty_Kernel_Time(i) = mean(var_Kernel_Time(coh == cohs(i)));
        avgRT_Kernel_Time(i) = mean(time_Kernel_Time(coh == cohs(i)));
    end
    if 1
        figure;
        plot(cohs, avgChoice_Kernel, 'r--'); hold on;
        plot(cohs, avgCertainty_Kernel, 'b--'); hold off;
        title('PPC with Correlation')
    end
    
    figure;
    plot(cohs, avgChoice_Kernel_Time, 'r--'); hold on;
    title(['Choice: PPC with Correlation and Time w/ Bound =' num2str(bound_Kernel_Time)])
    figure;
    plot(cohs, avgCertainty_Kernel_Time, 'b--'); hold off;
    title(['Variance: PPC with Correlation and Time w/ Bound =' num2str(bound_Kernel_Time)])
    figure;
    plot(cohs, avgRT_Kernel_Time, 'k--'); hold on; 
    title(['RT: PPC w/ Bound =' num2str(bound_Kernel_Time)])
end

%%
% Play around with the bound and see if it plays around with
% Confidence, Choice, and RT graphs -> It obviously should. (Break time
% more finely, and see if that changes anything -> instead of using 50ms
% use 1ms)
% - It be cool if you can back and did the two pool method Wei Ji Ma does
% for two differenct populations. It would test how different set of
% populations work with PPC

%% some sanity checks


% posterior as a function of coherence
% I dont think that is the posterior, i thought that was the likelihod
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
