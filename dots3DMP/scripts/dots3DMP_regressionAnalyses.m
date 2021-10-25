% function dots3DMP_regressionAnalyses(data)

% Regression models

% SJ 06/2020 started
% SJ 10/2021 modified significantly

% Analysis of behavioral data

% 1. improvement in accuracy when combined condition is presented
% 2. effect of heading strength and RT on high bet probability/SEP
% (logistic for PDW, multiple linear reg for SEP)
% 3. effect of heading strength and RT on probability correct (logistic)
% 4. improvement in accuracy when high bet is chosen
tstatfun = @(beta,se) 1-tcdf(beta./se,length(beta)-1);


mods = unique(data.modality);
cohs = unique(data.coherence);
deltas = unique(data.delta);

%% 1. improvement in accuracy with combined condition

J = data.delta==0;
I = data.modality(J)==3;
hdg = abs(data.heading(J));

D = [ones(size(hdg)), I, hdg.*I hdg, data.correct(J)];
% D = [ones(size(data.heading)), data.PDW, (data.heading).*data.PDW (data.heading), data.choice-1];
[beta,llik,pred,se] = logistfit_se(D);
p = tstatfun(beta,se);


%% 2. effect of heading strength on confidence

J = true(size(data.heading));
% J = data.modality==mods(2) & data.coherence==cohs(1);

if conftask == 1
    
elseif conftask == 2
    % P(High Bet) = [1 + exp(-b0 + b1*hdg + b2*RT) ] ^ -1
    
    D = [ones(size(data.heading(J))), abs(data.heading(J)), data.PDW(J)];
%     D = [ones(size(data.heading)), data.RT, abs(data.heading), data.PDW];
    [beta,llik,pred,se] = logistfit_se(D);
    
    p = tstatfun(beta,se);
    
end

% p(3) is <<0.05, and coefficient is positive, i.e. stronger heading
% increases probability of choosing high bet

%% 3. effect of heading strength on P(correct)

% P(correct) = [1 + exp(-b0 + b1*hdg + b2*RT) ] ^ -1
    
D = [ones(size(data.heading)), abs(data.heading), data.correct];
[beta,llik,pred,se] = logistfit_se(D);
p = tstatfun(beta,se);

% beta +ve and p significant for both RT and heading --> longer RT does
% result in increased p(correct) (on average), stronger heading increases
% p(correct) as well

%% 3b. effect of confidence on P(correct)

% is this valid given logistfit_se structure?? or must heading always be
% the last indep variable in D

D = [ones(size(data.heading)), abs(data.PDW), data.correct];
[beta,llik,pred,se] = logistfit_se(D);
p = tstatfun(beta,se);

%% 4. examine whether higher confidence is associated with increased accuracy

D = [ones(size(data.heading)), data.PDW, abs(data.heading).*data.PDW abs(data.heading), data.correct];
% D = [ones(size(data.heading)), data.PDW, (data.heading).*data.PDW (data.heading), data.choice-1];
[beta,llik,pred,se] = logistfit_se(D);
p = tstatfun(beta,se);

% if we use 
% p(2:4) all significant, indicating p(high bet) improves accuracy across
% all headings

%% 5.are confidence judgements associated with different reaction times

% separately for each mod and coh!


p = nan(length(mods),length(cohs),3);
for m=1:length(mods)
    for c=1:length(cohs)      
        J = data.modality==mods(m) & data.coherence == cohs(c) & data.delta == 0;
        
        if sum(J)==0, continue, end
        [p(m,c,:),t,stats,terms] = anovan(data.RT(J),[abs(data.heading(J)) data.PDW(J)],'model','interaction','display','off');
    end
end


%%

%{



% Fit regression models to choice, RT, conf data
% Choice data is going to be fitted by logistic regression (MLE method
% under binomial assumptions i.e. trials are a list of Bernoulli)




cohs   = unique(data.coherence);
hdgs   = unique(data.heading);
deltas = unique(data.delta);
mods   = unique(data.modality); 


tstatfun = @(beta,se) 1-tcdf(beta./se,length(beta)-1);


%%
% 1. quantify change in sensitivity associated with combined condition

I = data.modality==3;

% Pright = 1 + exp(-(b0 + b1*I + b2*I*hdg + b3*hdg))
D = [ones(size(data.heading)), I, I.*data.heading, data.heading data.choice-1];
%                                                               ^ subtract 1 to make choices 0 or 1                                                                

[beta,llik,pred,se] = logistfit_se(D);

% null hypothesis is that combined condition has no effect, i.e. b2 = 0
% tstat = (m-mu)/se  (mu=0 so m-mu is just beta)
p = tstatfun(beta,se);

% reject null if p(3) < 0.05 --> combined condition interacts with hdg to
% affect choice, and beta coefficient is positive, so combined condition significantly increases sensitivity  
% compare results when indicator variable refers to ves or vis only!

%% RT and conf (c.f. Kiani 2009/2014 analyses)

% 2. fit relationship between heading, RT and conf, separately for correct and
% error trials

goodRT = true(size(data.RT));

% goodRT = data.RT~=mode(data.RT); % kluge, in hand-tuned sim there is a risk that a
% lot of RTs will max out at stimulus time, so that too many weak stim (and
% error) trials have the same fixed RT - exclude these for now

% normalize, esp if there are different ranges across sessions/subjects
%RT = data.RT/max(data.RT); 
RT = data.RT;

% excluding zero heading trials, if any
data.correct = (data.choice==2 & data.heading>0) | (data.choice==1 & data.heading<0);
cor = data.correct & goodRT; 
err = ~data.correct & goodRT;

% correct trials
S = data.conf(cor);
% D = [ones(size(data.heading(cor))), abs(data.heading(cor)), RT(cor)];
% [Bcor,Bint_cor,~,~,statscor] = regress(S,D);

% use fitlm instead of regress (Bcor should match Estimate from fitlm)
tbl = table(abs(data.heading(cor)),RT(cor),S,'VariableNames',{'Hdg','RT','Conf'});
lmCor = fitlm(tbl,'Conf~Hdg+RT');

% repeat for error trials
S = data.conf(err);
% D = [ones(size(data.heading(err))), abs(data.heading(err)), RT(err)];
% [Berr,Bint_err] = regress(S,D);

tbl = table(abs(data.heading(err)),RT(err),S,'VariableNames',{'Hdg','RT','Conf'});
lmErr = fitlm(tbl,'Conf~Hdg+RT');

% subjects are more confident for faster decisions, both correct and
% incorrect


%% 3. does slope of regression change for correct and error trials?

% null hypothesis, b5 = 0
%I = abs(data.heading)<=2.5; % use only low hdg trials
I = abs(data.heading)>0.01; % use all trials
I = I & goodRT;
S = data.conf(I); 

% D = [ones(size(data.heading(I))), abs(data.heading(I)), RT(I),...
%         err(I), err(I).*abs(data.heading(I)), err(I).*RT(I)];
% [Bcorerr,Bint_corerr] = regress(S,D);

% or use fitlm
tbl = table(abs(data.heading(I)),RT(I),err(I),S,'VariableNames',{'Hdg','RT','Corr','Conf'});
lmCorErr = fitlm(tbl,'Conf~Hdg+RT+Corr+Corr*Hdg+Corr*RT');
 
% no - RT:Corr interaction coefficient is not different from 0

%% 4. relationship between hdg and certainty for errors
% I = abs(data.heading)<=6; % use only low hdg trials (not enough after RT exclusion!)
S = data.conf(err & I); 

% D = [ones(size(data.heading(err & I))), abs(data.heading(err & I))];
% [Bhdgerr,Bint_hdgerr] = regress(S,D);

tbl = table(abs(data.heading(err&I)),S,'VariableNames',{'Hdg','Conf'});
lmHdgErr = fitlm(tbl,'Conf~Hdg');

% significant +ve coefficient - on error trials, subjects are more certain
% for larger heading directions


%% 5. effect of RT and hdg on P(correct)

% Pcor= 1 / (1 + exp(-(b0 + b1*RT + b2*hdg)))
D = [ones(size(data.heading(goodRT))), RT(goodRT), abs(data.heading(goodRT)), data.correct(goodRT)];
[beta,llik,pred,se] = logistfit_se(D);

% null hypothesis is that RT has no effect on Pcor , i.e. b1 = 0
p = 1-tcdf(beta./se,length(beta)-1);

%% does modality (comb vs ves) influence the effect of RT on p(correct)?

% this is wrong...

temp = data;
removethese = data.modality==2;
for F = 1:length(fnames)
    eval(['temp.' fnames{F} '(removethese) = [];']);
end
RTtemp = RT;
RTtemp(removethese) = [];

D = [ones(size(temp.heading)), RTtemp, temp.modality==3,...
    (temp.modality==3).*abs(temp.heading), temp.correct];
% D = [ones(size(temp.heading)), RTtemp, temp.modality==3, abs(temp.heading),...
%     (temp.modality==3).*abs(temp.heading), temp.correct];
[beta,llik,pred,se] = logistfit_se(D);

% null hypothesis is that modality has no influence on the relationship between RT and Pcor , i.e. b4 = 0
p = 1-tcdf(beta./se,length(beta)-1);



%%%
%%
% Behavioral Analyses

% SJ 05/2021
% monkey PDW task

clear;clc

% load 'clean' data struct
load('lucio_20210315-20210707_clean.mat')

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



%}
