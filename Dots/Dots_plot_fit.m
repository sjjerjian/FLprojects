function Dots_plot_fit(parsedData,parsedFit,cohs,conftask,RTtask)

% somewhat unwieldy amalgam of parseData and Dots_plot, to show model 
% fits versus data in Dots experiment

% Spawned from dots3DMP_plots_fit_byCoh, but much simpler than that file,
% because it gets parsedData and parsedFit as inputs, and can call the
% existing Dots_plot for the data

% CF 01/2022

if nargin<3, conftask=0; end
if nargin<4, RTtask = 0; end


%% plot the data

wFit = 1;
Dots_plot(parsedData,cohs,conftask,RTtask,wFit)


%% plot the fit

nplots = 1+RTtask+double(conftask>0);

% all trials
figure(101);
subplot(nplots,1,1); hold on;
plot(cohs, parsedFit.pRight,'k-'); hold on;
if RTtask
    subplot(nplots,1,2); hold on;
    plot(cohs, parsedFit.RT, 'b-');
    if conftask
        subplot(nplots,1,3); hold on;
        plot(cohs, parsedFit.pHigh, 'r-');
    end    
else 
    if conftask
        subplot(nplots,1,2); hold on;
        plot(cohs, parsedFit.RT, 'b-');
    end
end


% separate choice+RT by high/low bet, and p(high) by corr/err
figure(102);
subplot(nplots,1,1); hold on;
plot(cohs, parsedFit.pRightHigh, 'k-');
plot(cohs, parsedFit.pRightLow, 'k--');
if RTtask
    subplot(nplots,1,2); hold on;
    plot(cohs, parsedFit.RTmeanHigh, 'b-');
    plot(cohs, parsedFit.RTmeanLow, 'b--');
    if conftask
        subplot(nplots,1,3); hold on;
        plot(cohs, parsedFit.pHighCorr, 'r-');
        plot(cohs, parsedFit.pHighErr, 'r--');
    end
else
    if conftask
        subplot(nplots,1,2); hold on;
        plot(cohs, parsedFit.pHighCorr, 'r-');
        plot(cohs, parsedFit.pHighErr, 'r--');
    end
end





