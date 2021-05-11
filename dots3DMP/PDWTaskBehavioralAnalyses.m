% Behavioral Analyses

% SJ 05/2021
% monkey PDW task

clear;clc

% load 'clean' data struct
load('lucio_20210401-20210505_proc.mat')

%%
% 1. quantify change in sensitivity associated with combined condition

I = data.modality==3;

% Pright = 1 + exp(-(b0 + b1*I + b2*I*hdg + b3*hdg))
D = [ones(size(data.heading)), I, I.*data.heading, data.heading data.choice-1];
%                                                               ^ subtract 1 to make choices 0 or 1                                                                
[beta,llik,pred,se] = logistfit_se(D);

% tstat = (m-mu)/se  (mu=0 so m-mu is just beta)
pChoice = 1-tcdf(beta./se,length(beta)-1);

% null hypothesis is that combined condition has no effect, i.e. b2 = 0
% reject null as p(3) (b2) << 0.05 --> combined condition interacts with hdg to
% affect choice, and beta coefficient is positive, so combined condition significantly increases sensitivity  
% (compare results when indicator variable I refers to ves or vis only)

%% 2. effect of heading on...
% probability of high bet

goodRT = data.RT<=1.5;
[n,edges,RTbin] = histcounts(data.RT(goodRT),10);

D = [ones(size(data.heading)), abs(data.heading), data.PDW];
[beta,llik,pred,se] = logistfit_se(D);
pHigh = 1-tcdf(beta./se,length(beta)-1);

% and probability correct
D = [ones(size(data.heading)), abs(data.heading), data.correct];
[beta,llik,pred,se] = logistfit_se(D);
pCor = 1-tcdf(beta./se,length(beta)-1);

%% 3. improvement in accuracy on high bets

D = [ones(size(data.heading)), abs(data.heading), data.PDW, data.correct];
[beta,llik,pred,se] = logistfit_se(D);
pHighAcc = 1-tcdf(beta./se,length(beta)-1);


%% 4. effect of heading and ves/comb modality on RT  

I = data.modality==3;
D = [ones(size(data.heading)), abs(data.heading), I, abs(data.heading).*I];

X = data.modality~=2;% & abs(data.heading)<=2.5;
D = D(X,:); % ignore visual
RT = data.RT(X);

[B,Berr,~,~,stats] = regress(RT,D);

tbl = table(D(:,2),D(:,3),RT,'VariableNames',{'Hdg','Modality','RT'});
lmHdgMod = fitlm(tbl,'RT~Hdg+Modality+Hdg*Modality');




% TO DO 

% a way to incorporate RT into logistic functions

% design plots of probability high bet and probability correct as a
% function of RT quantiles (a la Kiani & Shadlen 2009 Figure 1)

