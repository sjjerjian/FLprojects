function [wvesEmp,wvesPred] = dots3DMP_wgts_thres_NN(muPMF,sigmaPMF,cohs,deltas)

D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

for c = 1:length(cohs)      % m c d
    wvesPred(c) = (1/sigmaPMF(1,1,D)^2) / ((1/sigmaPMF(1,1,D)^2) + (1/sigmaPMF(2,c,D)^2));
                     % m c d
    actual(1) = (muPMF(3,c,1)-muPMF(3,c,2)+(deltas(1)/2)) / deltas(1);
    actual(2) = (muPMF(3,c,3)-muPMF(3,c,2)+(deltas(3)/2)) / deltas(3);    
    wvesEmp(c) = mean(actual);
end
