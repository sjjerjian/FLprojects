%% Tutorial on Binomial/Multinomial, and Beta/Dirichlet distributions
% SJJ, 07/2020

% The Dirichlet distribution is related to the multinomial distribution as the beta distribution is related to the binomial
% so let's start with a quick primer on the binomial-beta distributions and their relationship.

% The binomial distribution models the number of successes x in n trials.

% Imagine tossing a coin, and let x be the number of heads
% We can use the binomial distribution to model the probability of x heads
% The distribution has two parameters, n, the number of trials, and p, the probability of success
% x is thus a random variable, p (and n) is a parameter

% generate n random numbers from binomial distribution
% i.e. we flip a coin n times, count how many heads, then repeat N times
% R is the distribution of heads that we get
N = 10000;
n = 100;
p = 0.4; % for an unbiased coin, p=0.5
R = binornd(n,p,N,1);

% expected value and variance
% mu  = n*p; var = n*p*(1-p);
[mu,var] = binostat(n,p);
sd  = sqrt(var); % standard deviation

% from this simulation, we can empirically estimate the probability of
% getting x heads, just by summing
x = 40;
Px_hist = sum(R==x)/N

% if we assume the underlying distribution for our data is binomial 
% (here we know that it is, because we set it to be binomial in the first place!), 
% we can use the actual probability density function. 
Px_pdf = binopdf(x,n,p)

% Px_hist and Px_pdf should get closer and closer together as N increases

% plot it
figure; subplot(211); hold on
hh=histogram(R);
plot([mu mu],[0 N],'k','linew',2)
plot([mu-sd mu-sd],[0 N],'k--','linew',1.25) % mean-1SD
plot([mu+sd mu+sd],[0 N],'k--','linew',1.25) % mean+1SD

h_pdf = plot(binopdf(1:n,n,p)*N,'r','linew',2); % theoretical pdf
axis([0 n 0 N/8])

xlabel('number of successes (heads)');
ylabel('runs');
legend([hh,h_pdf],'histogram','pdf')
title(sprintf('Binomial distribution, p = %.1f, n = %d',p,n));

%% the beta distribution 

% the beta distribution provides a distribution over p, i.e. probability of success is now the random variable
% there are two parameters, alpha and beta, controlling shape of the
% distribution (see https://en.wikipedia.org/wiki/Beta_distribution to see
% what happens as you change the parameters)

% we can imagine a is the probability of success, b is the probability of failure
% so if we observe 40% success in a binomial process, what is the probability that the underlying
% p of the binomial distribution is 0.4? 

a = n*p; b = n*(1-p); 
Rbeta = betarnd(a,b,N,1);

subplot(212); hold on
histogram(Rbeta);
xlabel('P');
ylabel('count');
title('Beta distribution, a = 40, b = 60')

plot(0:0.01:1,betapdf(0:0.01:1,a,b)*n,'r','linew',2) % theoretical pdf, scaled
axis([0 1 0 N/8])

% the beta distribution is a distribution over a probability distribution (it tells us the probability
% of each value of p in the original binomial)

% Using the histogram data from the binomial distribution, we can estimate
% the parameters of the beta distribution
phat = betafit(R/n)
% ^ pretty close!

%% The multinomial and Dirichlet distributions

% these simply generalize the above to multiple outcomes (>2)

% multinomial distribution, where we can have k categories
n  = 10000;
p  = [0.45 0.35 0.2]; % p is now a k-length vector 
Rm = mnrnd(n,p);

% Rm is the number of trials on which we will will observe each outcome.
% Like the binomial, as n gets larger, Rm/n should get closer to p.

% the dirichlet distribution is a distribution over each value in p, just
% the way the beta distribution was a distribution over p from the binomial

% now we can use the gamma distribution to plot a distribution around each value in p 
% (or similarly, a histogram of multiple gamrnd samples)
clear Rm_gampdf
x = linspace(1,n,n/10);
for c=1:length(Rm)
    Rm_gampdf(:,c) = gampdf(x,Rm(c),1);
end

% Rm gives us the 'shape' of the distribution, and b=1 is the scale
% so we can use this to take a sample from the Dirichlet distribution with
% that shape, for each value in p
Rm_Dirichlet = gamrnd(Rm,1);
Rm_Dirichlet = Rm_Dirichlet./sum(Rm_Dirichlet)*n; % normalize so that sum is still n

% therefore, we can also get a similar result by sampling 
% just run mnrnd many times, and then empirically estimate the distribution associated with each p
N = 1000;
Rm_sample = mnrnd(n,p,N);

figure('color','w'); hold on; 
for c=1:length(Rm)
    histogram(Rm_sample(:,c),'facecolor','w','normalization','pdf');
end
hg = plot(x,Rm_gampdf,'linew',2); 
legend(hg,strcat('p = ', string(num2cell(p))),'fontsize',12);

