function [X, LL_final, data, fit] = fitDDM_wConfidence_simple(D,options,guess,fixed)

if options.feedback==2
    options.fh = figure;
    hold on;   
    xlabel('Call number' );
    ylabel('-LL');
end
    
% parameter bounds for fitting
%      [k   B  theta alpha]
X_lo = [0.2 10 0.8 0 ];
X_hi = [0.7 40  3 0.3];


global call_num; call_num=1;

if all(fixed)
    X = guess;
    [LL_final, data] = err_fcn_combined_wConf(X,D,options);
else

    switch options.fitMethod
        case 'fms'
            fitOptions = optimset('Display', 'final', 'MaxFunEvals', 500*sum(fixed==0), 'MaxIter', ... 
                500*sum(fixed==0), 'TolX', 1e-3, 'TolFun', 1e-2, 'UseParallel', 'Always');
            [X, LL_final] = fminsearch(@(x) err_fcn_combined_wConf(x,D,options), guess, fitOptions);

        case 'global'
            % GlobalSearch from Global Optimization Toolbox
            fitOptions = optimoptions(@fmincon,'Display','iter',...
                'Algorithm','interior-point',...
                'FinDiffType','central',...
                'UseParallel','always');
            problem = createOptimProblem('fmincon','x0',guess(fixed==0),'objective',...
                @(x) err_fcn_combined_wConf(x,D,options),'lb',X_lo(fixed==0),...
                'ub',X_hi(fixed==0),'options',fitOptions);
            gs = GlobalSearch;
            [X,LL_final,flagg,outptg,manyminsg] = run(gs,problem);      
        case 'multi'
            % MultiStart from Global Optimization Toolbox
            fitOptions = optimoptions(@fmincon,'Display','iter',...
                'Algorithm','interior-point',...
                'FinDiffType','central',...
                'UseParallel','always');
            problem = createOptimProblem('fmincon','x0',guess(fixed==0),'objective',...
                @(x) err_fcn_combined_wConf(x,D,options),'lb',X_lo(fixed==0),...
                'ub',X_hi(fixed==0),'options',fitOptions);
            ms = MultiStart('UseParallel','always');
            [X,LL_final,flagm,outptm,manyminsm] = run(ms,problem,200);    
        case 'pattern'
            % Pattern Search from Global Optimization Toolbox
            fitOptions = psoptimset(@patternsearch);
            fitOptions = psoptimset(fitOptions,'CompletePoll','on','Vectorized','off',...
                'UseParallel','always');
            [X,LL_final,exitflag,output] = patternsearch(@(x) err_fcn_combined_wConf(x,D,options),...
                guess(fixed==0),[],[],[],[],X_lo(fixed==0),X_hi(fixed==0),[],fitOptions);
    end

end


% to generate smooth curves, call it again with interpolated cohs and sampled durs
ndursamples = 1000;
nstrsteps = 33; % should be odd so there's a zero
g_str = linspace(min(data.strength),max(data.strength),nstrsteps);
rand_dur = randsample(data.dur, ndursamples, 'true');

fit = struct;
fit.strength = repmat(g_str', ndursamples, 1);
for j = 1:ndursamples
    fit.dur((j-1)*nstrsteps+1 : (j-1)*nstrsteps+nstrsteps) = rand_dur(j);
end
fit.dur = fit.dur';
fit.choice = ones(size(fit.strength));
fit.conf = ones(size(fit.strength));

options.plot = 0; options.feedback=0;
[~,fit] = err_fcn_combined_wConf(X,fit,options);


end



