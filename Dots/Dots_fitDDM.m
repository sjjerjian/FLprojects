function [X, LL_final, fit, parsedFit] = Dots_fitDDM(guess,fixed,data,options)

% master function for Dots DDM fitting (3DMP remains separate)
%   data:  data struct, requires at minimum a variable for choice and one for signed coherence
%   options: options struct; includes the objective (error) function, i.e. which model to fit
%   guess: initial guess of params
%   fixed: which params are fixed versus vary in the optimization

disp('fitting...');

if options.feedback==2
    options.fh = figure;
    hold on;   
    xlabel('Call number' );
    ylabel('-LL');
end
    
% upper and lower bounds on the parameters
LB = guess/4;
UB = guess*4;

global call_num; call_num=1;

if all(fixed)
    X = guess;
else

    switch options.fitMethod
        case 'fms'
            % fminsearch, pretty standard; main alternatives are fmincon or fminunc
            fitOptions = optimset('Display', 'final', 'MaxFunEvals', 500*sum(fixed==0), 'MaxIter', ... 
                500*sum(fixed==0), 'TolX', 1e-3, 'TolFun', 1e-2, 'UseParallel', 'Always');
            [X,~] = fminsearch(@(x) options.errfcn(x,guess,fixed,data,options), guess, fitOptions);

        case 'global'
            % GlobalSearch from Global Optimization Toolbox
            fitOptions = optimoptions(@fmincon,'Display','iter',...
                'Algorithm','interior-point',...
                'FinDiffType','central',...
                'UseParallel','always');
            problem = createOptimProblem('fmincon','x0',guess(fixed==0),'objective',...
                @(x) options.errfcn(x,guess,fixed,data,options),'lb',LB(fixed==0),...
                'ub',UB(fixed==0),'options',fitOptions);
            gs = GlobalSearch;
            [X,~,~,~,~] = run(gs,problem);      
        case 'multi'
            % MultiStart from Global Optimization Toolbox
            fitOptions = optimoptions(@fmincon,'Display','iter',...
                'Algorithm','interior-point',...
                'FinDiffType','central',...
                'UseParallel','always');
            problem = createOptimProblem('fmincon','x0',guess(fixed==0),'objective',...
                @(x) options.errfcn(x,guess,fixed,data,options),'lb',LB(fixed==0),...
                'ub',UB(fixed==0),'options',fitOptions);
            ms = MultiStart('UseParallel','always');
            [X,~,~,~,~] = run(ms,problem,200);    
        case 'pattern'
            % Pattern Search from Global Optimization Toolbox
            fitOptions = psoptimset(@patternsearch);
            fitOptions = psoptimset(fitOptions,'CompletePoll','on','Vectorized','off',...
                'UseParallel','always');
            [X,~,~,~] = patternsearch(@(x) options.errfcn(x,guess,fixed,data,options),...
                guess(fixed==0),[],[],[],[],LB(fixed==0),UB(fixed==0),[],fitOptions);
        case 'bads'
            % Bayesian adaptive direct search; work in progress: see dots3DMP_fitDDM
            
            % "plausible" lower/upper bounds
%             PLB = guess/2;
%             PUB = guess*2;

            error('BADS code not ready');
    end


end

% run err func again at the fitted/fixed params to generate a final
% error value and model-generated data points (expected Pright etc)
options.ploterr = 0;
[LL_final, fit, parsedFit] = options.errfcn(X,guess,fixed,data,options);

disp('done.');

end



