% Behavioral Analyses

% SJ 05/2021
% monkey PDW task

clear;clc

% load 'clean' data struct
load('lucio_20210401-20210505_clean.mat')

%% 1. effect of heading on...
% a. probability of high bet

% goodRT = data.RT<=1.5;
% [n,edges,RTbin] = histcounts(data.RT(goodRT),10);

D = [ones(size(data.heading)), abs(data.heading), data.PDW];
[beta,llik,pred,se] = logistfit_se(D);
pHigh = 1-tcdf(beta./se,length(beta)-1);

% b. probability correct
D = [ones(size(data.heading)), abs(data.heading), data.correct];
[beta,llik,pred,se] = logistfit_se(D);
pCor = 1-tcdf(beta./se,length(beta)-1);

%% 2. improvement in accuracy on high bets (sort of the same above, but both together)

D = [ones(size(data.heading)),  data.PDW, abs(data.heading), data.correct];
[beta,llik,pred,se] = logistfit_se(D);
pHighAcc = 1-tcdf(beta./se,length(beta)-1);


%% 4. effect of heading and ves/comb modality on RT  
% not sure about this

I = data.modality==3;
D = [ones(size(data.heading)), I, abs(data.heading).*I, abs(data.heading)];

X = data.modality~=2;% & abs(data.heading)<=2.5;
D = D(X,:); % ignore visual
RT = data.RT(X);

[B,Berr,~,~,stats] = regress(RT,D);

tbl = table(D(:,2),D(:,3),RT,'VariableNames',{'Hdg','Modality','RT'});
lmHdgMod = fitlm(tbl,'RT~Hdg+Modality+Hdg*Modality');


%% 5. effect of delta on confidence (low headings only)

%K = abs(data.heading)<5;
D = [ones(size(data.heading)), abs(data.delta) data.PDW];
%D = D(K,:);
[beta,llik,pred,se] = logistfit_se(D);
pDelta = 1-tcdf(beta./se,length(beta)-1);

% TO DO 

% a way to incorporate RT into logistic functions

% design plots of probability high bet and probability correct as a
% function of RT quantiles (a la Kiani & Shadlen 2009 Figure 1)

