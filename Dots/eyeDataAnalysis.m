% eyeDataAnalysisWrapper.m

clear
load Hanzo_behav_eyedata


%%

ppd = 37.9084; % pixels per degree (eye data are in pixel units, but targs are in degrees)
sc = [959.50 539.50]; % screen center (pixel coordinates)

result = cell(length(numfile),1);

for n=1:length(numfile)
    
    nTr = length(numfile{n}.eyeXYs);
    result{n} = nan(nTr,1);
    
    % plot/animate some individual trials:
    for Tr = 1:nTr
        
        if numfile{n}.twoTarget(Tr)==0 || numfile{n}.twotargconfidence(Tr)==0
            continue
        end
        
        % convert to degrees, with [0,0] as the center
        eyeX = (numfile{n}.eyeXYs{Tr}(1,:) - sc(1)) / ppd;
        eyeY = (numfile{n}.eyeXYs{Tr}(2,:) - sc(2))*-1 / ppd; % y pixels are inverted, because origin is upper-left
        % remove the NaNs
        eyeX(isnan(eyeX))=[];
        eyeY(isnan(eyeY))=[];
        if length(eyeX)~=length(eyeY); error('X and Y are different lengths'); end

        % one style: plot x and y as separate time series
        figure(11); clf;
        taxis = 1/120:1/120:length(eyeX)/120;
        plot(taxis,eyeX,'b',taxis,eyeY,'r');
        ylabel('deg'); xlabel('time (s)'); legend('X','Y');
        
        % another style: plot eye 'cursor' as x,y coordinate
        figure(12); set(gcf,'Color',[1 1 1],'Position',[800 900 670 380]); clf;
%         plot([-25 25],[0 0],'k--',[0 0],[-14 14],'k--'); hold on;
        plot(0,0,'k.','MarkerSize',12); hold on;
        plot(numfile{n}.locationLeftLowTarget(1),numfile{n}.locationLeftLowTarget(2),'b.','MarkerSize',16);
        plot(numfile{n}.locationRightLowTarget(1),numfile{n}.locationRightLowTarget(2),'b.','MarkerSize',16);
        plot(numfile{n}.locationLeftHighTarget(1),numfile{n}.locationLeftHighTarget(2),'g.','MarkerSize',16);
        plot(numfile{n}.locationRightHighTarget(1),numfile{n}.locationRightHighTarget(2),'g.','MarkerSize',16);
        circle(numfile{n}.aperture{Tr}(1),numfile{n}.aperture{Tr}(2),numfile{n}.aperture{Tr}(3)/2);            
        set(gca,'xlim',[-25 25],'ylim',[-14 14]); title(['trial ' num2str(Tr)]);
        h = nan(1,length(5:firstNaN(Tr)-1));
        pause(0.5);
        for t = 5:firstNaN(Tr)-1
            h(t) = plot(eyeX(t-4:t),eyeY(t-4:t),'rx'); % can change this to show more time points as desired
            pause(0.02); % controls the speed of the animation
            if t>5; set(h(5:t-1),'visible','off'); end
        end
        
        
        
        % Next task is to come up with a way to detect a saccade (based on
        % filtering or taking the derivative of the eye position data), use
        % it to identify the  subset of trials where the monkey made more
        % than one saccade (a putative change-of-mind or -of-confidence),
        % and then use the animations to convince yourself that most of 
        % these trials look like you'd expect. 
        
        % Once fairly convinced, confirm some expected patterns:
        % 1) are there more changes-of-mind when coherence is low?
        % what about when confidence is low, for a given coherence?
        % 2) do changes of mind improve accuracy? (more often changing
        % from incorrect to correct),
        % 3) see Resulaj et al. or van den Berg et al. for others
        
        % OR, start with calculating velocity (vigor) and plot it as a
        % function of conf or accuracy. Then might have more confidence
        % that your velocity metric (which has to be calculated from the
        % proper time window) is correct and reliable for detecting changes
        % of mind etc.
        

        
        %        result{n}(Tr) = whatever; %placeholder
        
        
    end

 end
