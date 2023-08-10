function ax = PTE_plotCueW(vesW, vesW_avg, clr)
    numColumns = size(vesW{1}, 2);
    numRows = size(vesW{1}, 1);
%     numel(vesW)
%     length(vesW)
%     size(vesW,2)

    figure;
    hold on;
    for c = 1:length(vesW)
        % Plot data from vesW
        for col = 1:numColumns
            x = col * ones(numRows, 1) + c/10;
            y = vesW{c}(:,col);
            plot(x, y, 'o', 'MarkerSize', 7, 'MarkerEdgeColor', clr{c});
        end

        % Plot data from vesW_avg
        for col = 1:numColumns
            x = col * ones(1, 1) + c/10 ;
            y = vesW_avg{c}(:,col);
            plot(x, y, 'o', 'MarkerSize', 6, 'MarkerFaceColor', clr{c}, 'MarkerEdgeColor', clr{c});
        end
    end
    hold off;
    ax = gca;
    ax.XTick = [1 + ((c-1)/2)/10, 2 + ((c-1)/2)/10];
    ax.XTickLabels = {'High coh trials','Low coh trials'};
    title('Previous trial coherence and modality affect cue weighting x2'); ylabel('Vestibular cue weight');
    ax.XLim = [.75,2.45];
end