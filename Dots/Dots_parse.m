cohs = unique(data.scoh);
dirs = unique(data.direction);

n = nan(length(cohs),1);

n1 = n; % all
n2 = n; % high
n3 = n; % low
nRT1 = n; % all
nRT2 = n; % high
nRT3 = n; % low
nPDW1 = n; % all
nPDW2 = n; % corr
nPDW3 = n; % err

pRight = n1;
pRightHigh = n2;
pRightLow = n3;

RTmean = nRT1;
RTse = nRT1;
RTmeanHigh = nRT2;
RTseHigh = nRT2;
RTmeanLow = nRT3;
RTseLow = nRT3;

pHigh = nPDW1;
pHighCorr = nPDW2;
pHighErr = nPDW3;

% for logistic regression (or whatever)
xVals = linspace(cohs(1),cohs(end),100);

for c = 1:length(cohs)
    % choice
    J = data.scoh==cohs(c);
    n1(c) = sum(J);
    pRight(c) = sum(J & data.choice==2) / n1(c); % 1 is left, 2 is right
    
    JJ = data.scoh==cohs(c) & data.PDW==1;
    n2(c) = sum(JJ);
    pRightHigh(c) = sum(JJ & data.choice==2) / n2(c); 
    
    JJJ = data.scoh==cohs(c) & data.PDW==0;
    n3(c) = sum(JJJ);
    pRightLow(c) = sum(JJJ & data.choice==2) / n3(c); 
    
    % RT
    K = ~isnan(data.RT); % need to start tagging non-RT trials as NaN!
    nRT1(c) = sum(J & K);
    RTmean(c) = mean(data.RT(J & K));
    RTse(c) = std(data.RT(J & K))/sqrt(nRT1(c));
    
    nRT2(c) = sum(JJ & K);
    RTmeanHigh(c) = mean(data.RT(JJ & K));
    RTseHigh(c) = std(data.RT(JJ & K))/sqrt(nRT2(c));

    nRT3(c) = sum(JJJ & K);
    RTmeanLow(c) = mean(data.RT(JJJ & K));
    RTseLow(c) = std(data.RT(JJJ & K))/sqrt(nRT3(c));
    
    % conf
    L = ~isnan(data.PDW);
    nPDW1(c) = sum(J & L);
    pHigh(c) = sum(J & L & data.PDW==1) / nPDW1(c); % 1 is high-bet

    LL = ~isnan(data.PDW) & data.correct==1;
    nPDW2(c) = sum(J & LL);
    pHighCorr(c) = sum(J & LL & data.PDW==1) / nPDW2(c);

    LLL = ~isnan(data.PDW) & data.correct==0;
    nPDW3(c) = sum(J & LLL);
    pHighErr(c) = sum(J & LLL & data.PDW==1) / nPDW3(c);
end

pRightSE = sqrt( (pRight.*(1-pRight)) ./ n1 );
pRightSEhigh = sqrt( (pRightHigh.*(1-pRightHigh)) ./ n2 );
pRightSElow = sqrt( (pRightLow.*(1-pRightLow)) ./ n3 );

pHighSE = sqrt( (pHigh.*(1-pHigh)) ./ nPDW1 );
pHighSEcorr = sqrt( (pHighCorr.*(1-pHighCorr)) ./ nPDW2 );
pHighSEerr = sqrt( (pHighErr.*(1-pHighErr)) ./ nPDW3 );


% fit logistic regression
X = data.scoh;
y = data.choice==2; % 2 is right
[B1, ~, stats1] = glmfit(X, y, 'binomial');

yVals1 = glmval(B1,xVals,'logit');

I = data.PDW==1;
X = data.scoh(I);
y = data.choice(I)==2; % 2 is right
[B2, ~, stats2] = glmfit(X, y, 'binomial');
yVals2 = glmval(B2,xVals,'logit');

I = data.PDW==0;
X = data.scoh(I);
y = data.choice(I)==2; % 2 is right
[B3, ~, stats3] = glmfit(X, y, 'binomial');
yVals3 = glmval(B3,xVals,'logit');