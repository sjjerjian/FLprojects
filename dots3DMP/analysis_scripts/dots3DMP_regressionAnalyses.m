% Regression models

% SJ 06/2020 started

% Fit regression models to choice, RT, conf data
% Choice data is going to be fitted by logistic regression (MLE method
% under binomial assumptions i.e. trials are a list of Bernoulli)

% note that confidence is treated as continuous variable here, i.e. this is
% for saccadic endpoint task only!

% simulations 
% folder = '/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/';
% file = 'RTdata_4-28-20.mat';

% file = 'simKiani09_20200506.mat';
% file = '2DAcc_simdata.mat';
% load([folder file],'data','cohs','deltas','hdgs','mods')

cohs   = unique(data.coherence);
hdgs   = unique(data.heading);
deltas = unique(data.delta);
mods   = unique(data.modality); 

% should already be removed
% fnames = fieldnames(data);
% removethese = data.oneConfTargTrial==1;
% for F = 1:length(fnames)
%     eval(['data.' fnames{F} '(removethese) = [];']);
% end

%%
% 1. quantify change in sensitivity associated with combined condition

I = data.modality==3;

% Pright = 1 + exp(-(b0 + b1*I + b2*I*hdg + b3*hdg))
D = [ones(size(data.heading)), I, I.*data.heading, data.heading data.choice-1];
%                                                               ^ subtract 1 to make choices 0 or 1                                                                
[beta,llik,pred,se] = logistfit_se(D);

% null hypothesis is that combined condition has no effect, i.e. b2 = 0
% tstat = (m-mu)/se  (mu=0 so m-mu is just beta)
p = 1-tcdf(beta./se,length(beta)-1);

% reject null as p(3) << 0.05 --> combined condition interacts with hdg to
% affect choice, and beta coefficient is positive, so combined condition significantly increases sensitivity  
% compare results when indicator variable refers to ves or vis only!

%% RT and conf (c.f. Kiani et al 2014 analyses)

% 2. fit relationship between heading, RT and conf, separately for correct and
% error trials

goodRT = true(size(data.RT));

% goodRT = data.RT<maxRT; % kluge, in hand-tuned sim there is a risk that a
% lot of RTs will max out at stimulus time, so that too many weak stim (and
% error) trials have the same fixed RT - exclude these for now

% normalize, esp if there are different ranges across sessions/subjects
RT = data.RT/max(data.RT); 
% RT = data.RT;

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


