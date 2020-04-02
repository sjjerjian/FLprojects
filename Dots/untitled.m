s
% % 
% % 
% % %     choice based on summed spikes from two pools (Shadlen et al. 1996)
% %     choiceP(t) = 2*(sum(r(prefpool)) > sum(r(nullpool))) - 1;
% %     
% %     % now try different ways to perturb R to see if we can degrade
% %     % performance of poisson decoder (but not simple pooling model)
% %     r = sum(squeeze(R(t,:,1:dur(t))),2) / (dur(t)/1000);
% % %     r(prefpool) = r(prefpool) + dur(t)/1000*50;
% % %     r(nullpool) = r(nullpool) + dur(t)/1000*50;
% %     r = r + 30*(sign(randn(size(prefpool))).*prefpool)';
% %     r(r<=0)=0.01;
% %     h = log(tuning{c}/10);
% %     L = exp(transpose(h)*(r/10));
% %     Lnorm = L/sum(L);
% %     lik2(t,:) = Lnorm;
% %     choiceL2(t) = sign(log(L(1)/L(181)));
% %     choiceP2(t) = 2*(sum(r(prefpool)) > sum(r(nullpool))) - 1;





% MultiStart from Global Optimization Toolbox
fitOptions = optimoptions(@fmincon,'Display','iter',...
    'Algorithm','interior-point',...
    'FinDiffType','central',...
    'UseParallel','always');

    beta_lo = beta_guess*0.1;
    beta_hi = beta_guess*10;
    problem = createOptimProblem('fmincon','x0',beta_guess,'objective',...
        @(j) err_fcn(j,dirAxis,Lnorm'),'lb',beta_lo,'ub',beta_hi,'options',fitOptions);
    ms = MultiStart('UseParallel','always');
    [X,fval,flagm,outptm,manyminsm] = run(ms,problem,200);    
