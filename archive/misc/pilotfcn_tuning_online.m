function pilotfcn_tuning_online(data)

subtractBaseline = 0;
frWindow = [100 ceil(max(data.dur)*1000)];
bkgWindow = [-50 20];
PSTHoffset = 50; % offset for PSTH matrix, e.g. 50 means plot starts at -50 ms

convKernel = fspecial('average',[1 50]);

loc_ = [data.ap_x, data.ap_y];
size_ = data.ap_diam/10;    %in degrees
speed_ = data.speed/10;     %in degrees/sec
dir_ = data.dir;

fr = nan(size(data.spikes,1),1);
bkg = nan(size(data.spikes,1),1);
trial_raster = nan(size(data.spikes,1),frWindow(2)+50+PSTHoffset);
for n = 1:size(data.spikes,1)
    spikesAlignedOnset = round((data.spikes{n}-data.dotsOn(n))*1000);
    I = spikesAlignedOnset >= frWindow(1) & spikesAlignedOnset <= frWindow(2);
    fr(n) = sum(I) * 1000/(frWindow(2)-frWindow(1));
    I = spikesAlignedOnset >= bkgWindow(1) & spikesAlignedOnset <= bkgWindow(2);
    bkg(n) = sum(I) * 1000/(bkgWindow(2)-bkgWindow(1));
    spikesForPSTH = spikesAlignedOnset + PSTHoffset;
    trial_raster(n,:) = 0;
    trial_raster(n,spikesForPSTH(spikesForPSTH>0)) = 1;
    trial_raster(n,round(data.dur(n)*1000)+PSTHoffset+1:end) = NaN;
        % ^can't keep spikes after dots off because trial ends and no more
        % are counted, leading to artifically reduced rates at end of PSTH
end
   %subtract bakground rate from fr
if subtractBaseline
    fr = fr - bkg;
end

%********** End, calculate the firing rates *********

%********** Begin, calculate & draw tuning curves *********

% first calculate FR and PSTH on catch trials, then remove them
frCatch = nanmean(fr(isnan(dir_)));
frCatch_se = std(fr(isnan(dir_)))/sqrt(sum(isnan(dir_)));
psthCatch = nanmean(trial_raster(isnan(dir_),:))*1e3;
psthCatch = smoothRaster(psthCatch, convKernel);

fr(isnan(dir_)) = [];
trial_raster(isnan(dir_),:) = [];
loc_(isnan(dir_)) = [];
size_(isnan(dir_)) = [];
speed_(isnan(dir_)) = [];
dir_(isnan(dir_)) = [];

dir_set = unique(dir_);
if length(dir_set)>1
        %find the mean and se response for each direction
    [r, r_se] = calc_mean(fr, dir_, dir_set);
    dir_set = mod(dir_set, 360);
    [dir_set, I] = sort(dir_set);
    r = r(I);
    r_se = r_se(I);
        %find the tuning parameters 
    [R,I] = max(r);
    modelParam = find_prefDir(fr, dir_, [dir_set(I),R,1,min(r)], [0 0 0 0]);        
    if modelParam.final(2) > 0
        prefDir = modelParam.final(1);
    else
        prefDir = mod(modelParam.final(1)+180, 360);
    end 
        % draw the tuning curve 
    g_dir = 1:360;
    g_r = modelParam.final(2)*exp(modelParam.final(3)*cos((modelParam.final(1)-g_dir)*pi/180)) + modelParam.final(4);
    sfigure(5); clf; subplot(2,1,1);
    h = polar([dir_set; dir_set(1)]*pi/180, [r; r(1)]); hold on;
    set(h, 'LineWidth', 2); % changeAxesFontSize(gca, 14, 14);
    [x, y] = pol2cart(prefDir*pi/180,max(r));
    line([0 x], [0 y], 'Color', 'r');
    title(sprintf('Preferred direction = %3.1f',prefDir));
    h = polar(0:0.2:2*pi,ones(size(0:0.2:2*pi))*frCatch,'r-');
    set(h, 'LineWidth', 2);
    
    subplot(2,1,2);
    errorbar(dir_set, r, r_se, 'o-', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    hold on;
%     plot(g_dir, g_r, 'r');    
    plot(dir_set, ones(size(dir_set))*frCatch,'r-','LineWidth', 2);
    box off;
    set(gca, 'XLim', [min(dir_set)-5 max(dir_set)+5], 'TickDir', 'out');
    xlabel('Direction (degree)');
    ylabel('Firing rate (spikes/s)');
%     set(gca, 'XLim', [0 360], 'XTick', 0:45:360, 'XTickLabel', makeTickLabel(0:45:360,90), ...
%          'YLim',[-55 150], 'YTick', -50:50:150, 'YTickLabel', makeTickLabel(-50:50:150,50), 'TickDir', 'out');
%     changeAxesFontSize(gca, 16, 16);
    
        %show PSTHs for different directions
    psth = calc_mean(trial_raster, dir_, dir_set)*1e3;
    psth = smoothRaster(psth, convKernel);
    YLim = [min(psth(:))*0.9 max(psth(:))*1.1];  
    t = -PSTHoffset:size(trial_raster,2)-PSTHoffset-1;
    sfigure(6); clf;
    for i = 1 : min(length(dir_set),8)
        subplot(3,4,i); 
        plot(t, psth(i,:), 'k');
        set(gca, 'XLim', [t(1) t(end)], 'YLim', YLim);
        box off;
        title(sprintf('%3.1f deg',dir_set(i)));
        if i == 5
            xlabel('Time (ms)');
            ylabel('Spikes/s');
        end
    end
    subplot(3,4,9); title('catch');
    plot(t, psthCatch, 'r');
    set(gca, 'XLim', [t(1) t(end)], 'YLim', YLim);
    box off; xlabel('Time (ms)'); ylabel('Spikes/s');
    
end


size_set = unique(size_);
if length(size_set)>1
        %draw the size tuning curve
    [r, r_se] = calc_mean(fr, size_, size_set);
    sfigure(5);
    set(gcf, 'Color', 'w', 'Position', [50 100 340 230], 'PaperPositionMode', 'auto');
    errorbar(size_set, r, r_se, 'o-', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    box off;
    set(gca, 'XLim', [min(size_set) max(size_set)], 'TickDir', 'out');
    xlabel('Size (degree)');
    ylabel('Firing rate (spikes/s)');
    
        %show PSTHs for different RF size
    psth = calc_mean(trial_raster, size_, size_set)*1e3;
    psth = smoothRaster(psth, convKernel);
    YLim = [min(psth(:))*0.9 max(psth(:))*1.1];  
    t = psth_alignto.start_offset:psth_alignto.end_offset;
    figure;
    set(gcf, 'Color', 'w', 'Position', [400 100 880 480], 'PaperPositionMode', 'auto');
    for i = 1 : min(length(size_set),8)
        subplot(2,4,i); 
        hold on;
        plot(t, psth(i,:), 'k');
        set(gca, 'XLim', [t(1) t(end)], 'YLim', YLim);
        box off;
        title(sprintf('%3.1f deg',size_set(i)));
        if i == 5
            xlabel('Time (ms)');
            ylabel('Spikes/s');
        end
    end
end

speed_set = unique(speed_);
if length(speed_set)>1
        %draw the speed tuning curve
    [r, r_se] = calc_mean(fr, speed_, speed_set);
            %calc tuning parameters
    [R,I] = max(r); 
    options = optimset ( 'Display' , 'final' , 'MaxFunEvals' , 500*4 , 'MaxIter' , 500*4 );
    param = fminsearch ( @(x) tuning_speed_err(x,fr,speed_) , [speed_set(I),R,0.1,min(r)] , options );    
    prefSpeed = param(1);
            %draw the tuning curve     
    g_speed = 1:0.1:30;
    g_r = param(2)*exp(-param(3)*(log(param(1))-log(g_speed)).^2)+param(4);
    figure;
    set(gcf, 'Color', 'w', 'Position', [50 100 340 230], 'PaperPositionMode', 'auto');
    hold on;
    errorbar(log(speed_set), r, r_se, 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    plot(log(g_speed), g_r);
    box off;
    set(gca, 'XLim', [0.9*min(log(speed_set)) 1.1*max(log(speed_set))], 'XTick', log([2 3 4 6 8 12 18]), 'XTickLabel', [2 3 4 6 8 12 18], ...
        'TickDir', 'out');
    xlabel('Speed (degree/s)');
    ylabel('Firing rate (spikes/s)');
    title(sprintf('Preferred speed = %3.1f',prefSpeed));
    
        %show PSTHs for different speeds
    psth = calc_mean(trial_raster, speed_, speed_set)*1e3;
    psth = smoothRaster(psth, convKernel);
    YLim = [min(psth(:))*0.9 max(psth(:))*1.1];  
    t = psth_alignto.start_offset:psth_alignto.end_offset;
    figure;
    set(gcf, 'Color', 'w', 'Position', [400 100 880 480], 'PaperPositionMode', 'auto');
    for i = 1 : min(length(speed_set),8)
        subplot(2,4,i);        
        plot(t, psth(i,:), 'k');
        set(gca, 'XLim', [t(1) t(end)], 'YLim', YLim);
        box off;
        title(sprintf('%3.1f deg/sec',speed_set(i)));
        if i == 5
            xlabel('Time (ms)');
            ylabel('Spikes/s');
        end
    end
end

loc_set = unique(loc_,'rows');
if size(loc_set,1)>1
    r = calc_mean(fr, loc_, loc_set);
    figure;
    set(gcf, 'Color', 'w', 'Position', [100 100 440 330], 'PaperPositionMode', 'auto');
    hold on;
    for i = 1 : size(loc_set,1)
        C = [0.90 0.90 0.90] - (r(i)-min(r))/(max(r)-min(r))*0.90;
        rectangle('Position', [loc_set(i,1)-0.25 loc_set(i,2)-0.25 0.5 0.5], 'Curvature', [1 1], 'EdgeColor', C, 'FaceColor', C);
    end
    axis equal
    set(gca, 'TickDir', 'out');
     
        %fit a 2D gaussian 
    gauss2d = @(param,loc) param(3)*exp(-((loc(:,1)-param(1)).^2+(loc(:,2)-param(2)).^2)/param(4)^2)+param(5);
% % %     gauss2d = inline('param(3)*exp(-((loc(:,1)-param(1)).^2+(loc(:,2)-param(2)).^2)/param(4)^2)+param(5)', 'param', 'loc');
    guess_center = [loc_set(:,1)'*r, loc_set(:,2)'*r]/sum(r);
    guess_width = sqrt(guess_center*guess_center');
    guess_max = max(r);
    guess_min = min(r);
    param = lsqcurvefit(gauss2d, [guess_center, guess_max, guess_width, guess_min], loc_, fr);
    fprintf('center = [%3.1f %3.1f]\n width = %3.1f\n\n', param([1 2 4]));
    [g_x, g_y] = ndgrid(min(loc_set(:,1))-1:0.5:max(loc_set(:,1))+1, min(loc_set(:,2))-1:0.5:max(loc_set(:,2))+1);
    g_r = gauss2d(param, [g_x(:) g_y(:)]);
    g_r = reshape(g_r, size(g_x));
    figure;
    set(gcf, 'Color', 'w', 'Position', [100 100 440 330], 'PaperPositionMode', 'auto');
    hold on;
    pcolor(g_x, g_y, g_r);
    shading interp;
    title(sprintf('center = [%3.1f %3.1f]',param(1:2)));
    set(gca, 'TickDir', 'out');
end

%********** End, calculate & draw tuning curves *********



function err = tuning_speed_err(param, fr, speed)

prefSpeed = param(1);
maxResp = param(2);
widthParam = param(3);
baseResp = param(4);

mu = maxResp*exp(-widthParam*(log(prefSpeed)-log(speed)).^2)+baseResp;

err = sum((mu-fr).^2);


