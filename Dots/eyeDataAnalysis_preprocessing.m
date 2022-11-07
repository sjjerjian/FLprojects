% eyeDataAnalysisWrapper.m


% % temp, cleanup
% clear
% load Hanzo_behavioralData_AllOfIt
% for n=1:length(numfile)
%     numfile{n} = rmfield(numfile{n},'dotPos');
%     % numfile{n} = rmfield(numfile{n},'eyeFlipTimes');
%     try
%         numfile{n} = rmfield(numfile{n},'idDotPos');
%     catch
%     end
% end
% save Hanzo_behavioralData.mat numfile -v7.3


%%
clear
load Hanzo_behav_eyedata


%%
    
ppd = 37.9084; % pixels per degree (eye data are in pixel units, but targs are in degrees)
sc = [959.50 539.50]; % screen center (pixel coordinates)

for n=1:length(numfile)
    clear endEyeXYs RT
    % check whether eyeXYs turns to NaNs at the time of RT:
    for Tr = 1:length(numfile{n}.eyeXYs)
        try
            firstNaN(Tr) = find(isnan(numfile{n}.eyeXYs{Tr}(1,:)),1,'first');
            % are there any NaNs interspersed?
            if any(~isnan(numfile{n}.eyeXYs{Tr}(1,firstNaN(Tr):end)))
                disp('YES!'); % nope
            end
        catch
            % last nan is not the same as RT, because it's sampled at 120
            % Hz, and is offset by pre-fixation time (and penalty phase?)
            eyedataTimes = numfile{n}.eyeFlipTimes{Tr}(1,1:firstNaN(Tr))-numfile{n}.eyeFlipTimes{Tr}(1,1);
            eyeRT = find(abs(eyedataTimes-numfile{n}.pdwEntered(Tr)) == min(abs(eyedataTimes-numfile{n}.pdwEntered(Tr)))) - ...
                    find(abs(eyedataTimes-numfile{n}.timeFpEntered(Tr)) == min(abs(eyedataTimes-numfile{n}.timeFpEntered(Tr))));
            eyeRT = eyeRT*8.333;
            altRT = (numfile{n}.pdwEntered(Tr)-numfile{n}.timeFpEntered(Tr))*1000;
            % ^ these two are matched, but that's kinda trivial.
            % what we want is an eyeRT that corresponds to numfile{n}.RT,
            % but we can't get that currently because we don't have time of
            % dots onset or saccade onset (only FP onset and targ acquire).
        end
        RT(Tr) = round(numfile{n}.RT(Tr)*1000);
    end
    figure(n);plot(RT,firstNaN*1/120*1000,'x')
    xlim([0 11000]); ylim([0 11000]); axis square;
    xlabel('RT');ylabel('eyeNaN');
    % they don't match, but are correlated in three distinct clouds--odd..
    % Oh, I think I know why: sometimes the intertrial interval is included
    % and more importantly, the time prior to fixation is highly variable.
    
    % anyway, leave that alone for now,
    % let's plot/animate some individual trials:
    for Tr = 1:length(numfile{n}.eyeXYs)
        
        % convert to degrees, with [0,0] as the center
        eyeX = (numfile{n}.eyeXYs{Tr}(1,:) - sc(1)) / ppd;
        eyeY = (numfile{n}.eyeXYs{Tr}(2,:) - sc(2))*-1 / ppd; % y pixels are inverted, because origin is upper-left
        % remove the NaNs
        eyeX(isnan(eyeX))=[]; 
        eyeY(isnan(eyeY))=[];
        if length(eyeX)~=length(eyeY); error('X and Y are different lengths'); end

        % one style: plot x and y as separate time series
%         figure(11); clf;
%         plot(eyeX,'b'); hold on;
%         plot(eyeY,'r');
        
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
            pause(0.02);
            if t>5; set(h(5:t-1),'visible','off'); end
        end
        
        % Next task is to come up with a way to detect a saccade (based on
        % filtering or taking the derivative of the eye position data), use
        % it to identify the  subset of trials where the monkey made more
        % than one saccade (a putative change-of-mind or -of-confidence),
        % and then use the animations to convince yourself that most of 
        % these trials look like you'd expect. 
        
        % Once fairly convinced, confirm some expected patterns:
        % 1) are there more changes-of-mind when coherence is low,
        % 2) do a change of mind improve accuracy? (more often changing
        % from incorrect to correct),
        % 3) see Resulaj et al. or van den Berg et al. for others
        
    end

 end
