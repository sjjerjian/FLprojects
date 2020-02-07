clear all;
close all;
load hanzo_20190501-20190803

% change to signed coh
data.coherence(data.direction==180) = -data.coherence(data.direction==180);

removethese = data.duration>1 | isnan(data.PDW);
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

coh_set = {[-0.512 0.512], [-0.256 0.256], [-0.128 0.128], [-0.064 0.064], [-0.032 0.032], 0};  %how to group coherences

colorscheme = {'m-',   'b-',     'g-',        'k-',     'r-',    'c-'};

Q = 4; % one over the fraction of the data to use in the running mean
logscale = 1;

% Pcorr vs. dur, sep trace for each unsigned coh (nostim only)
figure(11); set(gcf, 'Color', [1 1 1], 'Position', [100 200 450 475], 'PaperPositionMode', 'auto'); clf;
figure(12); set(gcf, 'Color', [1 1 1], 'Position', [550 200 450 475], 'PaperPositionMode', 'auto'); clf;
for c = 1:length(coh_set)
    I = ismember(data.coherence,coh_set{c}); % & data.direction==180; 
    
    [dur1,Pcorr,se_dp1] = running_mean(data.duration(I), data.correct(I)==1, sum(I)/Q);
    [dur2,Phigh,se_dp2] = running_mean(data.duration(I), data.PDW(I)==1, sum(I)/Q);

    if logscale
        figure(11);
        startX = dur1(1); 
        h11(c) = semilogx(dur1(find(dur1>=startX,1,'first'):end), Pcorr(find(dur1>=startX,1,'first'):end), colorscheme{c}, 'LineWidth',2); hold on;
        figure(12);
        startX = dur2(1); 
        h12(c) = semilogx(dur2(find(dur2>=startX,1,'first'):end), Phigh(find(dur2>=startX,1,'first'):end), colorscheme{c}, 'LineWidth',2); hold on;
    else
        figure(11);
        h11(c) = plot(dur1, Pcorr, colorscheme{c}, 'LineWidth',2); hold on;
        figure(12);
        h12(c) = plot(dur2, Phigh, colorscheme{c}, 'LineWidth',2); hold on;
    end
end


figure(11);
xlabel('Viewing duration (s)');
ylabel('Probability correct');
changeAxesFontSize(gca, 24, 24);
if logscale
    set(gca, 'XLim', [0.09 1.2], 'XTick', [0.15 0.3 0.6 1.2], 'YLim', [0.4 1], 'YTick', 0:0.1:1, 'TickDir', 'out');
else
    set(gca, 'XLim', [0.07 1], 'XTick', 0:0.2:1, 'YLim',[0.4 1], 'YTick', 0:0.2:1, 'TickDir', 'out');
end
set(gca, 'box', 'off');
l = legend(h11,{'51' '25' '12' '6' '3' '0'},'orientation','horizontal','location','northoutside');
set(l,'fontsize',12,'box','off')

% 'YTickLabel', makeTickLabel(0:0.2:1,0.2),

figure(12);
xlabel('Viewing duration (s)');
ylabel('Probability high bet');
changeAxesFontSize(gca, 24, 24);
if logscale
    set(gca, 'XLim', [0.09 1.2], 'XTick', [0.15 0.3 0.6 1.2], 'YLim', [0.4 0.9], 'YTick', 0.4:0.1:0.9, 'TickDir', 'out');
else
    set(gca, 'XLim', [0.07 1], 'XTick', 0:0.2:1,'YLim', [0.45 0.9], 'YTick', 0.4:0.1:0.9, 'TickDir', 'out');
end
set(gca, 'box', 'off');
l = legend(h12,{'51' '25' '12' '6' '3' '0'},'orientation','horizontal','location','northoutside');
set(l,'fontsize',12,'box','off')
