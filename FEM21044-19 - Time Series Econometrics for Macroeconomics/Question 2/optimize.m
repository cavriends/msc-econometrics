function [optimized_params, optimized_LL] = optimize(likelihood, startingvalues, y, fmincon_bounds)
% This function optimizes the negative log-likelihood using three different
% optimization algorithms (fminunc, fmincon, fminsearch).

%% fminunc
clearvars options
options = optimset('fminunc');
options = optimset(options, 'MaxFunEvals', 1e+6);
options = optimset(options, 'MaxIter', 1e+6);
options = optimset(options, 'TolFun', 1e-6);
options = optimset(options, 'TolX', 1e-6);

[fminunc_params, fminunc_LL] = fminunc(likelihood,startingvalues,options,y);

%% fmincon
clearvars options
options = optimset('fmincon');
options = optimset(options, 'MaxFunEvals', 1e+6);
options = optimset(options, 'MaxIter', 1e+6);
options = optimset(options, 'TolFun', 1e-6);
options = optimset(options, 'TolX', 1e-6);

lb = fmincon_bounds(:,1);
ub = fmincon_bounds(:,2);

[fmincon_params,fmincon_LL] = fmincon(likelihood, startingvalues,[],[],[],[],lb,ub,[],options,y);

%% fminsearch
clearvars options
options = optimset('fminsearch');
options = optimset(options, 'MaxFunEvals', 1e+6);
options = optimset(options, 'MaxIter', 1e+6);
options = optimset(options, 'TolFun', 1e-6);
options = optimset(options, 'TolX', 1e-6);

[fminsearch_params, fminsearch_LL] = fminsearch(likelihood, startingvalues, options, y);

optimized_params = [fminunc_params,fmincon_params,fminsearch_params];
optimized_LL = [fminunc_LL,fmincon_LL,fminsearch_LL];

end