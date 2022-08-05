function parsedData = Dots_parseData(data,conftask,RTtask,RTCorrOnly)

% Following SJ, CF converted to function, for cleaner workspace


if nargin < 4, RTCorrOnly = 0; end
if nargin < 3, RTtask = 0; end
if nargin < 2, conftask = 0; end

if RTtask==0
    data.RT = nan(size(data.choice));
end
switch conftask
    case 0
        data.PDW = nan(size(data.choice));
        data.conf = nan(size(data.choice));
    case 1
        data.PDW = nan(size(data.choice));
    case 2
        data.conf= nan(size(data.choice));
end

if ~isfield(data,'PDW_preAlpha') % before alpha offset, for choice/RT splits
    data.PDW_preAlpha = data.PDW;
end

cohs = unique(data.scoh);
if cohs(end)==1; cohs(end)=[]; end % one or more files had a few stray 100% coh trials in there by mistake

n = nan(length(cohs),1);

n1 = n; % all
n2 = n; % high
n3 = n; % low
pRight = n;
pCorrect = n;
pRightHigh = n;
pRightLow = n;

nRT1 = n; % all
nRT2 = n; % high
nRT3 = n; % low
RTmean = n;
RTse = n;
RTmeanHigh = n;
RTseHigh = n;
RTmeanLow = n;
RTseLow = n;

nPDW1 = n; % all
nPDW2 = n; % corr
nPDW3 = n; % err
pHigh = n;
pHighCorr = n;
pHighErr = n;

nConf1 = n; % all
nConf2 = n; % corr
nConf3 = n; % err
confMean = n;
confSE = n;
confMeanCorr = n;
confSEcorr = n;
confMeanErr = n;
confSEerr = n;


% for logistic regression (or whatever)
xVals = linspace(cohs(1),cohs(end),100);

K = ~isnan(data.RT); % need to start tagging non-RT trials as NaN!
if RTCorrOnly; K = K & data.correct==1; end

for c = 1:length(cohs)
    % choice
    J = data.scoh==cohs(c);
    n1(c) = sum(J);
    pRight(c) = sum(J & data.choice==1) / n1(c); % 0 is left, 1 is right

    pCorrect(c) = sum(J & data.correct==1) / n1(c); 
    
    if conftask==1
        Jhi = data.scoh==cohs(c) & data.conf>=median(data.conf);
        Jlo = data.scoh==cohs(c) & data.conf<median(data.conf);
    elseif conftask==2
        Jhi = data.scoh==cohs(c) & data.PDW_preAlpha==1;
        Jlo = data.scoh==cohs(c) & data.PDW_preAlpha==0;
    end
    if conftask
        n2(c) = sum(Jhi);
        pRightHigh(c) = sum(Jhi & data.choice==1) / n2(c); 
        n3(c) = sum(Jlo);
        pRightLow(c) = sum(Jlo & data.choice==1) / n3(c);        
    end
    
    % RT
    nRT1(c) = sum(J & K);
    RTmean(c) = mean(data.RT(J & K));
    RTse(c) = std(data.RT(J & K))/sqrt(nRT1(c));
    
    if conftask
        nRT2(c) = sum(Jhi & K);
        RTmeanHigh(c) = mean(data.RT(Jhi & K));
        RTseHigh(c) = std(data.RT(Jhi & K))/sqrt(nRT2(c));

        nRT3(c) = sum(Jlo & K);
        RTmeanLow(c) = mean(data.RT(Jlo & K));
        RTseLow(c) = std(data.RT(Jlo & K))/sqrt(nRT3(c));
    end
    
    % pdw
    if conftask==2
        L = ~isnan(data.PDW);
        nPDW1(c) = sum(J & L);
        pHigh(c) = sum(J & L & data.PDW==1) / nPDW1(c); % 1 is high-bet

        Lcor = ~isnan(data.PDW) & data.correct==1;
        nPDW2(c) = sum(J & Lcor);
        pHighCorr(c) = sum(J & Lcor & data.PDW==1) / nPDW2(c);

        Lerr = ~isnan(data.PDW) & data.correct==0;
        nPDW3(c) = sum(J & Lerr);
        pHighErr(c) = sum(J & Lerr & data.PDW==1) / nPDW3(c);
    end    
    
    % conf
    if conftask==1
        M = ~isnan(data.conf);
        nConf1(c) = sum(J & M);
        confMean(c) = mean(data.conf(J & M));
        confSE(c) = std(data.conf(J & M))/sqrt(nConf1(c));

        MM = ~isnan(data.conf) & data.correct==1;
        nConf2(c) = sum(J & MM);
        confMeanCorr(c) = mean(data.conf(J & MM));
        confSEcorr(c) = std(data.conf(J & MM))/sqrt(nConf2(c));

        MMM = ~isnan(data.conf) & data.correct==0;
        nConf3(c) = sum(J & MMM);
        confMeanErr(c) = mean(data.conf(J & MMM));
        confSEerr(c) = std(data.conf(J & MMM))/sqrt(nConf3(c));
    end    
end

pRightSE = sqrt( (pRight.*(1-pRight)) ./ n1 );
pCorrectSE = sqrt( (pCorrect.*(1-pCorrect)) ./ n1 );
pRightSEhigh = sqrt( (pRightHigh.*(1-pRightHigh)) ./ n2 );
pRightSElow = sqrt( (pRightLow.*(1-pRightLow)) ./ n3 );

pHighSE = sqrt( (pHigh.*(1-pHigh)) ./ nPDW1 );
pHighSEcorr = sqrt( (pHighCorr.*(1-pHighCorr)) ./ nPDW2 );
pHighSEerr = sqrt( (pHighErr.*(1-pHighErr)) ./ nPDW3 );


% fit logistic regression
% all trials
X = data.scoh;
y = data.choice==1; % 1 is right
[B1, ~, stats1] = glmfit(X, y, 'binomial');
yVals1 = glmval(B1,xVals,'logit');
    % high conf/bet
if conftask==1
    I = data.conf>=median(data.conf);
elseif conftask==2
    I = data.PDW_preAlpha==1;
else
    I = false(size(data.PDW));
end
X = data.scoh(I);
y = data.choice(I)==1;
[B2, ~, stats2] = glmfit(X, y, 'binomial');
yVals2 = glmval(B2,xVals,'logit');
    % low conf/bet
if conftask==1
    I = data.conf<median(data.conf);
elseif conftask==2
    I = data.PDW_preAlpha==0;
else
    I = false(size(data.PDW));
end
X = data.scoh(I);
y = data.choice(I)==1;
[B3, ~, stats3] = glmfit(X, y, 'binomial');
yVals3 = glmval(B3,xVals,'logit');
     

parsedData = struct();
parsedData.n = n;
parsedData.pRight = pRight;
parsedData.pRightHigh = pRightHigh;
parsedData.pRightLow = pRightLow;
parsedData.pRightSE = pRightSE;
parsedData.pRightSEhigh = pRightSEhigh;
parsedData.pRightSElow = pRightSElow;
parsedData.pCorrect = pCorrect;
parsedData.pCorrectSE = pCorrectSE;
parsedData.xVals = xVals;
parsedData.yVals1 = yVals1;
parsedData.yVals2 = yVals2;
parsedData.yVals3 = yVals3;
parsedData.B1 = B1;
parsedData.B2 = B2;
parsedData.B3 = B3;
parsedData.stats1 = stats1;
parsedData.stats2 = stats2;
parsedData.stats3 = stats3;
parsedData.yVals1 = yVals1;
parsedData.yVals2 = yVals2;
parsedData.yVals3 = yVals3;

if conftask==1
    parsedData.confMean = confMean;
    parsedData.confSE = confSE;
    parsedData.confMeanCorr = confMeanCorr;
    parsedData.confSEcorr = confSEcorr;
    parsedData.confMeanErr = confMeanErr;
    parsedData.confSEerr = confSEerr;
elseif conftask==2
    parsedData.pHigh = pHigh;
    parsedData.pHighSE = pHighSE;
    parsedData.pHighCorr = pHighCorr;
    parsedData.pHighSEcorr = pHighSEcorr;
    parsedData.pHighErr = pHighErr;
    parsedData.pHighSEerr = pHighSEerr;
end

if RTtask
    parsedData.RTmean = RTmean;
    parsedData.RTse = RTse;
    parsedData.RTmeanHigh = RTmeanHigh;
    parsedData.RTseHigh = RTseHigh;
    parsedData.RTmeanLow = RTmeanLow;
    parsedData.RTseLow = RTseLow;
end





