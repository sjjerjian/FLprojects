function wves = dots3DMP_cueWeights(gfit,cohs,deltas,conftask,ifplot)

if nargin<5, ifplot = 0; end
if nargin<4, conftask = 0; end

D = length(deltas)+1;

for c = 1:length(cohs)      % m c d
    
    % choice-based
    wves.choice.pred(c) = (1/gfit.choice.sigma(1,1,D)^2) / ((1/gfit.choice.sigma(1,1,D)^2) + (1/gfit.choice.sigma(2,c,D)^2));
    
    if length(deltas)>1
                                % m c d
    actual(1) = (gfit.choice.mu(3,c,1)-gfit.choice.mu(3,c,2)+(deltas(1)/2)) / deltas(1);
    actual(2) = (gfit.choice.mu(3,c,3)-gfit.choice.mu(3,c,2)+(deltas(3)/2)) / deltas(3);    
    wves.choice.emp(c) = mean(actual);
    
    end

    % confidence-based
    if conftask
        
        wves.conf.pred(c) = (1/gfit.conf.sigma(1,1,D)^2) / ((1/gfit.conf.sigma(1,1,D)^2) + (1/gfit.conf.sigma(2,c,D)^2));
        
        if length(deltas)>1
        % m c d
        actual(1) = (gfit.conf.mu(3,c,1)-gfit.conf.mu(3,c,2)+(deltas(1)/2)) / deltas(1);
        actual(2) = (gfit.conf.mu(3,c,3)-gfit.conf.mu(3,c,2)+(deltas(3)/2)) / deltas(3);
        wves.conf.emp(c) = mean(actual);
        
        end
    end
    
    %{
    %%% RT-based is probably not very useful, don't bother for now
    % RT-based
    if RTtask
        wves.RT.pred(c) = (1/gfit.RT.sigma(1,1,D)^2) / ((1/gfit.RT.sigma(1,1,D)^2) + (1/gfit.RT.sigma(2,c,D)^2));
        
        if length(deltas)>1
        % m c d
        actual(1) = (gfit.RT.mu(3,c,1)-gfit.RT.mu(3,c,2)+(deltas(1)/2)) / deltas(1);
        actual(2) = (gfit.RT.mu(3,c,3)-gfit.RT.mu(3,c,2)+(deltas(3)/2)) / deltas(3);
        wves.RT.emp(c) = mean(actual);
        
        end
    end
    %}
end

if ifplot
    % NEEDS TIDYING
    figure(501); hold on;
    
    plot(wves.choice.emp,'marker','o','color','r');
    plot(wves.choice.pred,'marker','o','color','k');
    
    if conftask
    plot(wves.conf.emp,'marker','o','color','r','linestyle','--');
    plot(wves.conf.pred,'marker','o','color','k','linestyle','--');
    end
end
    