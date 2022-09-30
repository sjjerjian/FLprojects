% Simple tutorial on the likelihood function as it pertains to neural decoding.
% Here there is just one neuron and 3 stimuli, but can be extended as needed.
%   - C. Fetsch, Oct.2011

% HOW TO USE:
% 1) First time around, run step by step (i.e., "Run Section" one at a time)
% while reading the comments.
% 2) Next, try uncommenting line 21 or 22 and run whole script again to test
% a different tuning curve


%% 1) present an identical stimulus (theta=theta1) many times and count spikes

clear all; close all;
% The spike count, r, is a random variable with some distribution,
% known as the probability density function (PDF) P(r|theta=theta1)
% (more accurately the probability mass function, since r is discrete).
% This distribution is what we are estimating by presenting a stimulus many times.

%%% ignore these for the first pass; then try uncommenting each and re-running:
% tc = [20 30 7]; r_obs = 23; % peaked tuning w/ linear likelihood function (same as default)
tc = [10 20 30]; r_obs = 20; % linear tuning w/ peaked likelihood function (alternate)

addedNoise = 0; % another thing you can play with


theta = 1; % arbitrarily call theta1 "1" on the theta axis
for n = 1:100 % 100 trials (the loop is unnecessary, but serves to emphasize these are sequential measurements)
    % assume a Poisson process, but NOTE: this assumption normally does NOT 
    % come first, before the measurements.  It's only necessary in a simulation.
    if ~exist('tc','var')
        r(n,theta) = poissrnd(20) + addedNoise; % 20 is arbitrary, it's just the mean response that theta1 happens to elicit
    else
        r(n,theta) = poissrnd(tc(theta));
    end
end
% r(:,theta) = poissrnd(20,100,1);

figure(1); set(gcf, 'Color', [1 1 1], 'Position', [200 200 1000 700], 'PaperPositionMode', 'auto'); hold on
subplot(2,3,1);
hist(r(:,theta));
xlabel('spike count (r)');
ylabel('number of times each r was observed');
title('spike histogram for one particular stimulus, theta1')


%% 2) repeat for a couple other values of theta

% These are separate and independent; don't try to link them in any way (yet)

theta = 2;
for n = 1:100
    if ~exist('tc','var')
        r(n,theta) = poissrnd(30) + addedNoise; % 30 is arbitrary
    else
        r(n,theta) = poissrnd(tc(theta));
    end
end
theta = 3;
for n = 1:100
    if ~exist('tc','var')
        r(n,theta) = poissrnd(7) + addedNoise; % 7 is arbitrary
    else
        r(n,theta) = poissrnd(tc(theta)); 
    end
end

subplot(2,3,2);
x_r = 0:1:50; % bin centers for counting spikes (histogram)
h = histc(r,x_r);
plot(x_r,h); legend('theta=theta_1','theta=theta_2','theta=theta_3')
xlabel('spike count (r)');
ylabel('number of times each r was observed');
title('spike histograms for three values of theta');


%% 3) plot the mean response as a function of theta, i.e. the tuning curve

tc = mean(r); % tuning curve
subplot(2,3,3);
errorbar([1 2 3], mean(r), std(r));
set(gca, 'Xtick', [1 2 3]);
xlim([0.75 3.25]);
xlabel('theta');
ylabel('mean spike count (+/- SD)');
title('tuning curve');


%% 4) now ask, given one particular observed response, what is the likelihood of each theta?

% First we need to have run enough trials to get some idea of the shape of the PDFs
% (here, let's say we have reason to believe they fit a Poisson distribution).
% This is an important step! Technically you cannot* compute the likelihood
% unless you have a parametric model of P(r|theta).

% % % P(model | data ) ~= [P(data | model) * P(model)] / P(data

% *(However, you can use the 'histogram' method (see below),
%   but typically this requires too much data to be practical)

subplot(2,3,4); hold on;
x_r = 0:50; % bin centers for counting spikes (histogram)
for t = 1:3
    r_model(:,t) = poisspdf(0:50,tc(t));
end

xlabel('spike counts (r)');
ylabel('probability P(r|theta=theta_n)');
title('PDFs for the three values of theta (poisson PMFs)');
plot(x_r,r_model);

if ~exist('r_obs','var')
    r_obs = 23; % the observation (a spike count)
end

% Where does our observation lie on these PMFs?  Plot each as a colored dot:
plot(r_obs, r_model(x_r==r_obs,1), 'b.', 'MarkerSize', 20); % theta1
plot(r_obs, r_model(x_r==r_obs,2), 'g.', 'MarkerSize', 20); % theta2
plot(r_obs, r_model(x_r==r_obs,3), 'r.', 'MarkerSize', 20); % theta3
legend('theta=theta_1','theta=theta_2','theta=theta_3', 'r_o_b_s')


% histogram method
for t = 1:3
    r_hist(:,t) = poissrnd(tc(t),100000,1); % make many fake trials at each theta
end
h = histc(r_hist,x_r) / 100000; % change histogram to probability
plot(x_r,h,'--','LineWidth',2);



%% Misc notes:

% Assuming we have enough trials, these points of intersection correspond to the values of the
% _likelihood function_, the (relative) probability of each possible value of theta, given the observed r

% *In words, this sounds like it should be written P(theta|r=r_obs), but nobody writes it this way, 
% probably to avoid confusion with the posterior.
% The posterior density P(theta|r) is a probability distribution whereas the likelihood (P(r|theta)) is not.
% It is composed of values taken from several independent probability distributions (the colored dots), and so 
% it generally does not sum to 1.

% Of course, the posteror P(theta|r) also depends on the prior, P(theta),
% by Bayes' rule -- but we don't need to deal with either of those here.


%% 5) plot the likelihood function (P(r|theta) as a function of theta, fixing r)

% Since our noise model is Poisson, we know the equation for P(r|theta)
for t=1:3
    L(t) = exp(-tc(t))*tc(t)^(r_obs) / factorial(r_obs); % "analytic" likelihood (from Poisson PMF equation)
    Lh(t) = h(x_r==r_obs,t); % "histogram-based" likelihood, as in step 4 above
end

subplot(2,3,5); hold on;
plot(1,L(1),'b.','MarkerSize',20);
plot(2,L(2),'g.','MarkerSize',20);
plot(3,L(3),'r.','MarkerSize',20);
plot([1 2 3], Lh, 'ks--');  % Note the overlap between histogram and analytic methods (sanity check)
set(gca, 'Xtick', [1 2 3]);
xlabel('theta');
ylabel('Likelihood (P(r|theta) as a function of theta)');
title('likelihood function');
legend('analytic','','','histogram');

%********
% Note that the likelihood looks nothing like the tuning curve, even though
% they both are functions of theta. The mapping between them depends on the
% noise model (e.g. Poisson) and the observed response relative to the 
% tuning curve.
%********


