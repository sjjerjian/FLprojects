%% Entropy tutorial
% SJJ 07-2020

% entropy quantifies the level of 'information' associated with a random
% variable.

% Using the coin toss example, consider flipping a coin with some known
% probability of heads/tails.

% in extreme case, where p(heads) = 0 or p(heads) = 1, the entropy is zero,
% since we don't get any new information from tossing the coin more times
% (there is no uncertainty)

% on the other hand, if p(heads) = 0.5, each coin toss is maximally
% uncertain, and it is very difficult to predict the outcome of the next
% toss. Each toss gives us new information (1 'bit' to be precise).
% P(heads) in between 0.5 and 1 will give us levels of information varying
% between these two extremes.

% a simple coin toss example

p = linspace(0,1,100); % p(heads) = 0.01, 0.01, 0.02, ..., 1
p(2,:) = 1-p;
p(p==0) = eps; % log(0) is undefined, so just change zeros to v small number

% calculate h for all p element-wise (don't need to write a loop)
h = -sum(p.*log2(p)); % formula for entropy

% Note that if p is 0 or 1, p.*log2(p) will be 0

figure('position',[600 500 300 200],'color','w');
plot(p,h,'linew',2,'marker','o'); box off;
xlabel('P(heads)')
ylabel('H(X)');
set(gca,'fontsize',12);
title('Entropy of coin toss')

%% we can do the same thing for a multinomial distribution
% entropy is highest when p = 1/n where n is the number of options

% define some arbitrary probability distributions for n=3 multinomials
pm = [1 eps eps; 0.7 0.2 0.1; 0.5 0.3 0.2; 1/3 1/3 1/3
      eps 1 eps; 0.3 0.5 0.2; 0.1 0.3 0.5];
hm = -sum(pm.*log2(pm),2);

% [~,maxH] = max(hm);
% pm(maxH,:)

figure('position',[200 500 300 200],'color','w');

subplot(211); bar(pm); box off;
ylabel('P(X)');
xlim([0 length(hm)]+0.5)
set(gca,'xtick',1:length(hm))

subplot(212); plot(hm,'linew',2,'marker','o','linestyle','none'); box off;
xlabel('Distr')
ylabel('H(X)');
xlim([0 length(hm)]+0.5)
set(gca,'xtick',1:length(hm))

for i=1:length(hm)
    fprintf('P = [%.g %.g %.g], H = %.2f\n',pm(i,1),pm(i,2),pm(i,3),hm(i))
end

% Note how entropy is minimal when two of the outcomes have probability 0
% (eps) and one has probability 1 (equivalent to the 1/0 coin toss
% example), and entropy is maximal for a uniform distribution (all
% probabilities equal 1/3)

%% Repeat multinomial, but across many random examples

nDist = 1000;

pmR = zeros(nDist,3);
pmR(:,1) = rand(nDist,1);

for c=1:size(pmR,2)-1
    pmR(:,c) = randrange(0,1-sum(pmR,2),[nDist,1]);
end

pmR(:,end) = 1 - sum(pmR,2);
pmR(pmR==0) = eps;

hmR = -sum(pmR.*log2(pmR),2);

[inds,edges] = discretize(hmR,nDist);
cmap = cbrewer('seq','PuRd',nDist);

hm_cols = cmap(inds,:);

figure;
scatter(pmR(:,1),pmR(:,2),40,hm_cols,'filled')


