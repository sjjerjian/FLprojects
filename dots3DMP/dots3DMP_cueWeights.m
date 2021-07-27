function wves = dots3DMP_cueWeights(gfit,cohs,deltas)

D = length(deltas)+1;

for c = 1:length(cohs)      % m c d
    
    % choice-based
    wves.choice.pred(c) = (1/gfit.choice.sigma(1,1,D)^2) / ((1/gfit.choice.sigma(1,1,D)^2) + (1/gfit.choice.sigma(2,c,D)^2));
                     % m c d
    actual(1) = (gfit.choice.mu(3,c,1)-gfit.choice.mu(3,c,2)+(deltas(1)/2)) / deltas(1);
    actual(2) = (gfit.choice.mu(3,c,3)-gfit.choice.mu(3,c,2)+(deltas(3)/2)) / deltas(3);    
    wves.choice.emp(c) = mean(actual);

    % confidence-based
    wves.conf.pred(c) = (1/gfit.conf.sigma(1,1,D)^2) / ((1/gfit.conf.sigma(1,1,D)^2) + (1/gfit.conf.sigma(2,c,D)^2));
                     % m c d
    actual(1) = (gfit.conf.mu(3,c,1)-gfit.conf.mu(3,c,2)+(deltas(1)/2)) / deltas(1);
    actual(2) = (gfit.conf.mu(3,c,3)-gfit.conf.mu(3,c,2)+(deltas(3)/2)) / deltas(3);    
    wves.conf.emp(c) = mean(actual);
    
    % RT-based
    wves.RT.pred(c) = (1/gfit.RT.sigma(1,1,D)^2) / ((1/gfit.RT.sigma(1,1,D)^2) + (1/gfit.RT.sigma(2,c,D)^2));
                     % m c d
    actual(1) = (gfit.RT.mu(3,c,1)-gfit.RT.mu(3,c,2)+(deltas(1)/2)) / deltas(1);
    actual(2) = (gfit.RT.mu(3,c,3)-gfit.RT.mu(3,c,2)+(deltas(3)/2)) / deltas(3);    
    wves.RT.emp(c) = mean(actual);

end