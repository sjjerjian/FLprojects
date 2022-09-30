%% answer to question about gambler's ruin in the Markov chain section.
% Define the starting state.
s0 = [0 0 0 0 1 0 0 0 0]

% Create the transition matrix
p  = .5,  q = 1-p
P = diag(p*ones(8,1),1) + diag(q*ones(8,1),-1); 
P(1,:) = 0; P(1,1) = 1; P(end,:)=0; P(end,end)=1
P
% s0 * P^16

nmax = 200
n = [1:nmax]';

% Go through step by step and extract the probability of being in one of
% the absorbtion states
cdf = zeros(nmax,1);
for i = 1:nmax
    S = s0 * P^i;
    cdf(i) = S(1) + S(end);
end
% The pdf is the derivative of the CDF
pdf = diff([0;cdf]);
trapz(n,pdf)  % area should be 1
sum(pdf)      % same thing because dn = 1.

figure(1), clf
plot(cdf)
figure(2), clf
plot(n,pdf)

% Mean and var from moments.
meanStoppingTime = sum(n .* pdf)
varStoppingTime = sum(n.^2 .* pdf) - meanStoppingTime^2 

% Now let's do this analytically
P
addpath('/Users/mike/Documents/Matlab_mfiles/mfiles/')

!grep Kemeny /Users/mike/Documents/Matlab_mfiles/mfiles/*.m
help absorbingMarkovToCanonical

[C,I,O,R,Q] = absorbingMarkovToCanonical(P)
% Fundamental matrix
N = inv(eye(size(Q)) - Q)
% mean absorption time from 
middleIndex = 1 + (length(Q)-1)/2
% Mean stopping time from each transient state
tau = N * ones(length(N),1)
tau(middleIndex) % should be the same as meanStoppingTime
% Variance of stopping time from each transient state
tau2 = (2*N - eye(size(N)))*tau - tau.^2
tau2(middleIndex)  % should be same as varStoppingTime

% reflecting lower bound
p  = .5,  q = 1-p
nsteps = 15
P = diag(p*ones(nsteps-1,1),1) + diag(q*ones(nsteps-1,1),-1); 
P(1,:) = 0; P(1,1) = 1;     % absorbing bound
P(end,:)=0; P(end,end-1)=1    % reflecting bound
[C,I,O,R,Q] = absorbingMarkovToCanonical(P);
% Fundamental matrix
N = inv(eye(size(Q)) - Q);
% Mean stopping time from each transient state
tau = N * ones(length(N),1)
% Variance of stopping time from each transient state
tau2 = (2*N - eye(size(N)))*tau - tau.^2

cv = sqrt(tau2) ./ tau

% Brute force approach
s0 = zeros(1,nsteps);   % row vector
s0(2) = 1
s0 * P
nmax = 2000
n = [1:nmax]';

cdf = zeros(nmax,1);
for i = 1:nmax
    S = s0 * P^i;
    cdf(i) = S(1);
end
figure(1), clf
plot(cdf)

figure(2), clf
pdf = diff([0;cdf]);
plot(n,pdf)

sum(pdf)

meanStoppingTime = sum(n .* pdf)
varStoppingTime = sum(n.^2 .* pdf) - meanStoppingTime^2 
cv = sqrt(varStoppingTime) / meanStoppingTime


