function [fitvariables] = SingleGaussFit(x0,LB,UB,x,y)
%This is a bare-bones fitting function that fits data (x,y) with a double
%gaussian specified by x0,LB,UB

%% Options
options = optimset(@lsqcurvefit);
options.Display = 'off';
options.TolFun = 1e-5;
options.MaxFunEvals = 1e3;
options.MaxIter = 200;
options.FinDiffType = 'central';

%% Fitting Function
fit_gauss = @(p) p(1)*exp(-(x-p(2)).^2./(2*p(3)^2));
err_fit_gauss = @(v) fit_gauss(v)-y;

%% Fitting Process
fitvariables = lsqnonlin(err_fit_gauss,x0,LB,UB,options);
end

