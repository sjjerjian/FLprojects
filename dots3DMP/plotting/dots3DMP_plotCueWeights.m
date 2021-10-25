function dots3DMP_plotCueWeights(wves,wvesBoot,cohs,conftask)

figure(810); set(gcf,'Color',[1 1 1],'Position',[50 20 500 250],'PaperPositionMode','auto'); clf;
if nargin==5,subplot(121); end
hold on;

erPred = errorbar(cohs,wves.choice.pred,std(wvesBoot.choice.pred),'k-o','LineWidth',2,'MarkerSize',6,'MarkerFaceColor','k'); 
erEmp = errorbar(cohs,wves.choice.emp,std(wvesBoot.choice.emp),'r-o','LineWidth',2,'MarkerSize',6,'MarkerFaceColor','r');    

set(erPred,'Color','k');
set(erEmp,'Color','r');

text(cohs(1),0.5,'Choice-Pred','color','k')
text(cohs(1),0.45,'Choice-Emp','color','r')

if conftask
    erConf = errorbar(cohs,wves.conf.emp,std(wvesBoot.conf.emp),'b-o','LineWidth',2,'MarkerSize',6,'MarkerFaceColor','c');
    set(erConf,'Color','b')
    text(cohs(1),0.4,'Conf-Emp','color','b')
end


ylim([0.25 1.05]);
set(gca,'XTick',cohs,'Xlim',[cohs(1)-0.05 cohs(end)+0.05],'XTickLabel',cohs,'Ytick',0:0.25:1);
xlabel('Visual coherence'); ylabel('Vestibular weight');
changeAxesFontSize(gca,16,16); set(gca,'box','off'); tidyaxes;


% plot thresholds

% if nargin==5
%     cols = 'krb';
%     subplot(122);
%     hold on;
%     
%     D = 4;
% 
%     % empirical thresholds
%     choice_thres = squeeze(gfit.choice.sigma(:,:,D));
%     conf_thres   = squeeze(gfit.conf.sigma(:,:,D));
%     
%     % normalize to ves
%     choice_thres = choice_thres ./ choice_thres(1,:);
%     conf_thres = conf_thres ./ conf_thres(1,:);
%     
%     temp = choice_thres.^2;
% 
%     choice_thres_pred = squeeze( (temp(1,:).*temp(2,:)) ./ (temp(1,:) + temp(2,:)) );
%     
%     for m=1:size(choice_thres,1)
%         plot(cohs,choice_thres(m,:),cols(m))
%     end
% end
    
    
    
