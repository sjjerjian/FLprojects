function hf = dots3DMP_plot_wgts_bootstrap(wvesPred,wvesEmp,wvesEmpboot,wvesPredboot,sigmaPMFboot,gfit,cohs)

hf = figure('position',[600 500 350 450]);

nboots = size(wvesEmpboot,1);

% vestibular weights
subplot(211); hold on
wvesEmp_err = abs(repmat(wvesEmp,2,1) - prctile(wvesEmpboot,[5 95],1));
errorbar(cohs,wvesEmp,wvesEmp_err(1,:),wvesEmp_err(2,:),'linestyle','-','color','black',...
    'marker','o','markeredgecolor','black','markerfacecolor','black')

wvesPred_err = abs(repmat(wvesPred,2,1) - prctile(wvesPredboot,[5 95],1));
errorbar(cohs,wvesPred,wvesPred_err(1,:),wvesPred_err(2,:),'linestyle','--','color','black',...
    'marker','o','markeredgecolor','black','markerfacecolor','white')
legend('Observed','Optimal')
axis([10 65 0 1])
set(gca,'xtick',cohs);

% thresholds
nn = nan(nboots,length(cohs));
Tves = nn; Tvis = nn;
TcombEmp = nn;

% full data (for means), or just take mean of bootstrapped data?
Tves = gfit.sigmaPMF(1,:,4);
Tvis = gfit.sigmaPMF(2,:,4);
TcombEmp = gfit.sigmaPMF(3,:,2);
TcombPred = sqrt( (Tves.^2 .* Tvis.^2) ./ (Tves.^2 + Tvis.^2) );

% bootstrapped data
for n=1:nboots
    Tves_boot(n,:) = sigmaPMFboot{n}(1,:,4);
    Tvis_boot(n,:) = sigmaPMFboot{n}(2,:,4);
    TcombEmp_boot(n,:) = sigmaPMFboot{n}(3,:,2);
end
TcombPred_boot = sqrt( (Tves_boot.^2 .* Tvis_boot.^2) ./ (Tves_boot.^2 + Tvis_boot.^2) );

% normalize everything to vestibular
Tves0 = Tves;
Tves = Tves ./ Tves0;
Tvis = Tvis ./ Tves0;
TcombEmp = TcombEmp ./ Tves0;
TcombPred = TcombPred ./ Tves0;

Tves_boot = Tves_boot ./ Tves0;
Tvis_boot = Tvis_boot ./ Tves0;
TcombEmp_boot = TcombEmp_boot ./ Tves0;
TcombPred_boot = TcombPred_boot ./ Tves0;

% compute confidence intervals
% errorbar requires actual values, so calculate diff from mean
Tves_err     = abs(repmat(Tves,2,1) - prctile(Tves_boot,[5 95],1));
Tvis_err     = abs(repmat(Tvis,2,1) - prctile(Tvis_boot,[5 95],1));
TcombEmp_err = abs(repmat(TcombEmp,2,1) - prctile(TcombEmp_boot,[5 95],1));
TcombPred_err = abs(repmat(TcombPred,2,1) - prctile(TcombPred_boot,[5 95],1));

subplot(212); hold on
errorbar(cohs,Tves,Tves_err(1,:),Tves_err(2,:),'linestyle','-','color','black',...
    'marker','o','markeredgecolor','black','markerfacecolor','black')
errorbar(cohs,Tvis,Tvis_err(1,:),Tvis_err(2,:),'linestyle','-','color','red',...
    'marker','o','markeredgecolor','red','markerfacecolor','red')
errorbar(cohs,TcombEmp,TcombEmp_err(1,:),TcombEmp_err(2,:),'linestyle','-','color','cyan',...
    'marker','o','markeredgecolor','cyan','markerfacecolor','cyan')
errorbar(cohs,TcombPred,TcombPred_err(1,:),TcombPred_err(2,:),'linestyle','--','color','cyan',...
    'marker','o','markeredgecolor','cyan','markerfacecolor','white')
legend('Ves','Vis','Observed','Optimal')
axis([10 65 0 1.5])
set(gca,'xtick',cohs);