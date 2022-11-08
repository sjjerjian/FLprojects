% fiddle_infoSDT.m
xmin = -10
xmax = 10
dx = .001
x = [xmin:dx:xmax]';
% The d'=1 case
uR = .5;
uL = -.5;
yR = normpdf(x,uR,1);
yL = normpdf(x,uL,1);
plot(x,yR,x,yL)
trapz(x,yL)
trapz(x,-log2(yR).*yR)
trapz(x,-log2(yL).*yL)

pdf_d = normpdf(x,uR-uL,sqrt(2));
plot(x,pdf_d)
% The ideal observers answers "RIGHT" when diff > 0
% How often does the ideal observers answer correctly (assuming that the
% true answer is right)
pcor = 1 - normcdf(0,uR-uL,sqrt(2))
p = [.5 .5]'
Htot = sum(-log2(p).*p)
pcond = [pcor 1-pcor]
Hcond = sum(-log2(pcond).*pcond)
MInfo = Htot - Hcond

% Entropy of the difference distribution
% Before we show the stimulus, it has equal prob of being in either
% direction
dpre = 0.5 * (normpdf(x,uR-uL,sqrt(2)) + normpdf(x,uL-uR,sqrt(2)));
dpost = pdf_d;
% Entropy in response, unconditional on stimulus
Hpre = trapz(x,-log2(dpre).*dpre)
% Entropy in response, conditional on stimulus
Hpost = trapz(x,-log2(dpost).*dpost)
% Mutual information between response and stimulus
Hpre-Hpost

% I don't understand why we get more mutual information this way
% 
% I'll try recoding the variable as binary by setting the probability of
% negative vals to 0
bpre = (x>=0) .* dpre;
a = sum(dx * bpre);
bpre = bpre/a;
bpost = (x>=0) .* dpost;
a = sum(dx * bpost);
bpost = bpost/a;
Hpre = nansum(dx * -log2(bpre).*bpre)
Hpost = nansum(dx * -log2(bpost).*bpost)
Hpre - Hpost





