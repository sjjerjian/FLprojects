% Stats for neurobiology.
format short g

%% SE, Var and SD
%  What is the SD and what is the SE

mu = 50
sigma= 10
n = 100
x = normrnd(mu,sigma,n,1)
% x = poissrnd(mu,n,1)
% x = exprnd(mu,n,1)
% x = raylrnd(mu,n,1)
% x = unifrnd(10,mu,n,1)

figure(1), clf
[nn xx] = hist(x);
hb1 = bar(xx,nn/n,1)
axprefs(gca)
xlabel('x')
ylabel('Frequency of observation')


xmean = mean(x)
xs = std(x)

% What is the variance? 
samplevar = var(x)
popvar = var(x,1)
mean((x - xmean).^2)
mean(x.^2) - mean(x).^2
xs*xs
% How much variation is there about the mean?
cv = xs/xmean

% What is the standard error?
% Let's imagine that instead of 100 samples in one experiment, leading to one mean, 
% we have Nexp experiments, each with 100 samples, hence 50 estimates of the mean.
Nexp = 100
n = 20
x = normrnd(mu,sigma,n,Nexp);
size(x)
% The "population" mean should be mu
mean(x(:))
% The population stdev should be sigma
std(x(:))
% We have Nexp sample means
sampleMeans = mean(x)
% What is the mean of the sample means?
mean(mean(x))
% What is the stdev of the sample means?
std(mean(x))
figure(1),hold on
[nn xx] = hist(sampleMeans);
hb2 = bar(xx,nn/Nexp,1); set(hb2,'FaceColor','r','EdgeColor','r')

% Compare this to the standard error
stdErrs = std(x)/sqrt(n)

mean(std(x)/sqrt(n))



% Why look at the stdev to tell us about uncertainty?
% Suppose we turn up the amplifier.
x = normrnd(mu,sigma,n,1)
y = 10 * x;
figure(2),clf
hist(y)
ymean = mean(y)
ys = std(y)
% What is the variance?
ys2 = var(y)
mean((y-ymean).^2)
mean(y.^2) - mean(y).^2
% How much variation is there about the mean?
ys/ymean


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison of two means
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make two groups of (fake) data
% Set means and n
m1 = 6
m2 = 4
n1 = 50
n2 = 40
x1 = poissrnd(m1,n1,1)
x2 = poissrnd(m2,n2,1)
% look at the data
[nn1 v1] = hist(x1);
[nn2] = hist(x2,v1);
figure(2),clf, hold on
h = plot(v1,nn1/n1,'b',v1,nn2/n2,'r'), set(h,'LineWidth',3)
bar(v1,[nn1(:)/n1 nn2(:)/n2])
axprefs(gca)
xlabel('value')
ylabel('Relative frequency')

% what are the means and SE
[mean(x1) std(x1)/sqrt(n1)]
[mean(x2) std(x2)/sqrt(n2)]
% How big is the difference? Let's invent a statistic:
% the absolute value of the difference between the means
absMeanDiffObs = abs(mean(x1) - mean(x2))

% Now let's play a game called The Null Hypothesis.
% Under H0, x1 and x2 are random draws from a single population
% Let's assume we have the entire population. 
u = [x1; x2]   
% u stands for union. The first n1 elements of u are 
% our original x1. The remaining n2 elements are x2.

% Let's make a bunch of these unions, 
% only let's rearrange the order.
nsamples = 200
U = repmat(0,n1+n2,nsamples);
size(U)
for i = 1:nsamples
	U(:,i) = u(randperm(n1+n2));
end
% Now cut these permuted union distributions into two groups.
% Respect the original sizes, n1 and n2
Y1 = U(1:n1,:);
Y2 = U(n1+1:end, :);
size(Y1)
size(Y2)
% Let's look at the differences between the means of these fake
% distributions. Since we don't care about the sign, it's probably
% wise to look at the absolute values
absMeanDiffH0 = abs(mean(Y1) - mean(Y2))
[nn vabs] = hist(absMeanDiffH0,20);
figure(1),clf
bar(vabs, nn/nsamples)
axprefs(gca)
xlabel('Absolute value of difference between means')
ylabel('Relative frequency')
% What about the means in our data
absMeanDiffObs = abs(mean(x1) - mean(x2))
% How often did we see a value this big or bigger in our 
% simulation under H0?
n = sum(absMeanDiffH0 >= absMeanDiffObs)
% Or the probability is 
n/nsamples
% If this is 0, then we should say the p < 1/nsamples
(n+1) / nsamples

