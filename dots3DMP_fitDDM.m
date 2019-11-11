function [X, LL_final, data, fit] = dots3DMP_fitDDM(D,options,guess,fixed)


if options.feedback>0
    options.fh = figure;
    hold on;   
    xlabel('Call number' );
    ylabel('-LL');
end


% parameter bounds for fitting
%      [kves kvis B]
X_lo = [0.1 0.5 20 ];
X_hi = [0.7 2.0 100];


global call_num; call_num=1;

if all(fixed)
    X = guess;
else

switch options.fitMethod
    case 'fms'
        fitOptions = optimset('Display', 'final', 'MaxFunEvals', 500*sum(fixed==0), 'MaxIter', ... 
            500*sum(fixed==0), 'TolX', 1e-3, 'TolFun', 1e-2, 'UseParallel', 'Always');
        [X, fval] = fminsearch(@(x) dots3DMP_fitDDM_err_stimStyle(x,D,options), guess, fitOptions);

    case 'global'
        % GlobalSearch from Global Optimization Toolbox
        fitOptions = optimoptions(@fmincon,'Display','iter',...
            'Algorithm','interior-point',...
            'FinDiffType','central',...
            'UseParallel','always');
        problem = createOptimProblem('fmincon','x0',guess(fixed==0),'objective',...
            @(x) dots3DMP_fitDDM_err_stimStyle(x,D,options),'lb',X_lo(fixed==0),...
            'ub',X_hi(fixed==0),'options',fitOptions);
        gs = GlobalSearch;
        [X,fval,flagg,outptg,manyminsg] = run(gs,problem);      
    case 'multi'
        % MultiStart from Global Optimization Toolbox
        fitOptions = optimoptions(@fmincon,'Display','iter',...
            'Algorithm','interior-point',...
            'FinDiffType','central',...
            'UseParallel','always');
        problem = createOptimProblem('fmincon','x0',guess(fixed==0),'objective',...
            @(x) dots3DMP_fitDDM_err_stimStyle(x,D,options),'lb',X_lo(fixed==0),...
            'ub',X_hi(fixed==0),'options',fitOptions);
        ms = MultiStart('UseParallel','always');
        [X,fval,flagm,outptm,manyminsm] = run(ms,problem,200);    
    case 'pattern'
        % Pattern Search from Global Optimization Toolbox
        fitOptions = psoptimset(@patternsearch);
        fitOptions = psoptimset(fitOptions,'CompletePoll','on','Vectorized','off',...
            'UseParallel','always');
        [X,fval,exitflag,output] = patternsearch(@(x) dots3DMP_fitDDM_err_stimStyle(x,D,options),...
            guess(fixed==0),[],[],[],[],X_lo(fixed==0),X_hi(fixed==0),[],fitOptions);
end

end
    
options.feedback = 0;
[LL_final, data] = dots3DMP_fitDDM_err_simStyle(X,D,options);


% to generate smooth curves, run again with interpolated headings
nstrsteps = 33; % should be odd so there's a zero
g_str = linspace(min(data.heading),max(data.heading),nstrsteps);

fit = struct;
fit.heading = repmat(g_str',1,10);
fit.choice = ones(size(fit.heading));
fit.conf = ones(size(fit.heading));

options.plot = 0;
[~,fit] = dots3DMP_fitDDM_err(X,fit,options);


end



