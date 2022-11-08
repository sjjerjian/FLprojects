function Dots_plot_fit(parsedData,parsedFit,cohs,conftask,RTtask,forTalk)

% plot model fits on top of preexisting data plots, for Dots paradigm

% Spawned from dots3DMP_plots_fit_byCoh, but much simpler than that file,
% because it gets parsedData and parsedFit as inputs, and can call the
% existing Dots_plot for the data

% CF 01/2022

if nargin<3, conftask=0; end
if nargin<4, RTtask = 0; end


%% plot the data

wFit = 1;
Dots_plot(parsedData,cohs,conftask,RTtask,wFit,forTalk)


%% plot the fit

nplots = 1+RTtask+double(conftask>0);

% all trials
figure(101);
subplot(nplots,1,1); hold on;
plot(cohs, parsedFit.pRight,'k-'); hold on;
if RTtask
    subplot(nplots,1,2); hold on;
    plot(cohs, parsedFit.RTmean, 'b-');
    if conftask
        subplot(nplots,1,3); hold on;
        plot(cohs, parsedFit.pHigh, 'r-');
    end    
else 
    if conftask
        subplot(nplots,1,2); hold on;
        plot(cohs, parsedFit.RTmean, 'b-');
    end
end

% separate choice+RT by high/low bet, and p(high) by corr/err
figure(102);
subplot(nplots,1,1); hold on;
plot(cohs, parsedFit.pRightHigh, 'k-');
plot(cohs, parsedFit.pRightLow, 'k--');
legend('high conf','low conf'); % temp, doubtconf
if RTtask
    subplot(nplots,1,2); hold on;
    plot(cohs, parsedFit.RTmeanHigh, 'b-');
    plot(cohs, parsedFit.RTmeanLow, 'b--');
    legend('high conf','low conf'); % temp, doubtconf
    if conftask
        subplot(nplots,1,3); hold on;
        plot(cohs, parsedFit.pHighCorr, 'r-');
        plot(cohs, parsedFit.pHighErr, 'r--');
        legend('corr','err','Location','North'); % temp, doubtconf
    end
else
    if conftask
        subplot(nplots,1,2); hold on;
        plot(cohs, parsedFit.pHighCorr, 'r-');
        plot(cohs, parsedFit.pHighErr, 'r--');
        legend('corr','err','Location','North'); % temp, doubtconf
    end
end

%% 

if forTalk

% choice
figure(105); hold on;
plot(cohs, parsedFit.pRight,'k-', 'LineWidth', 2);
export_fig('dots_choice','-eps');
% choice split
figure(106); hold on;
plot(cohs, parsedFit.pRightHigh, 'k-', 'LineWidth', 2);
plot(cohs, parsedFit.pRightLow, 'k--', 'LineWidth', 2);
h = legend('high bet','low bet');
set(h,'FontSize',16); legend('boxoff');
export_fig('dots_choice_split','-eps');

if RTtask    
    % RT
    figure(107); hold on;
    plot(cohs, parsedFit.RTmean, 'b-', 'LineWidth', 2);
    export_fig('dots_RT','-eps');
    % RT split
    figure(108); hold on;
    plot(cohs, parsedFit.RTmeanHigh, 'b-', 'LineWidth', 2);
    plot(cohs, parsedFit.RTmeanLow, 'b--', 'LineWidth', 2);
%     h = legend('high conf','low conf'); % temp, doubtconf
    h = legend('high bet','low bet');
    set(h,'position',[0.6512 0.7609 0.3222 0.1516]);
    set(h,'FontSize',16); legend('boxoff');
    export_fig('dots_RT_split','-eps');
end

if conftask
    % PDW
    figure(109); hold on;
    plot(cohs, parsedFit.pHigh, 'r-', 'LineWidth', 2);
	export_fig('dots_conf','-eps');
    % PDW split
    figure(110); hold on;
    plot(cohs, parsedFit.pHighCorr, 'r-', 'LineWidth', 2);
    plot(cohs, parsedFit.pHighErr, 'r--', 'LineWidth', 2);
    h = legend('correct','error','Location','Southwest');
    set(h,'FontSize',16); legend('boxoff');
    export_fig('dots_conf_split','-eps');
end


end


