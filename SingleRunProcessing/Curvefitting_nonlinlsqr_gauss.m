function [centroids,SNR,centroid_error,fitvariables_out,residual] = Curvefitting_nonlinlsqr_gauss(onerow_data,...
    gaussfit_limits,row,noise)
%This function fits gaussian curves to a 2-D set of intensity values over a
%pixel range. This function utilizes the nonlinear least squares fitting
%method.

%Just doing 2 gaussian curves right now for G1 and G2.

%% Trim to only use the central part of the image for the curve fitting
rowcounter = 1:length(onerow_data);

%% Options
options = optimset(@lsqcurvefit);
options.Display = 'off';
options.TolFun = 1e-5;
options.MaxFunEvals = 1e3;
options.MaxIter = 200;
options.FinDiffType = 'central';

%% Bounds (seperate into g1 and g2)
notwidths = [1,2,4,5];
LB = [gaussfit_limits(1,notwidths)];
x0 = [gaussfit_limits(2,notwidths)];
UB = [gaussfit_limits(3,notwidths)];
widths = gaussfit_limits(2,[3,6])';

%% Fitting G1 and G2
    %together
    fit_gauss_double_gauss = @(p) (p(1)*exp(-(rowcounter-p(2)).^2 ./ (2*widths(1)^2)))+(p(3)*exp(-(rowcounter-p(4)).^2 ./ (2*widths(2)^2)));
    err_fit_gauss_g = @(v) fit_gauss_double_gauss(v)-onerow_data;
    %fit
    [fitvariables,~,residual,~,~,~,jacobian] = lsqnonlin(err_fit_gauss_g,x0,LB,UB,options);
    %subtract the fit from the rest of data processing

%% Confidence in fit
%centroids
ci_g = nlparci(fitvariables,residual,'jacobian',jacobian);   %  95% confidence intervals for the fit coefficients
centroid_error = [(ci_g(2,2)-ci_g(2,1))/2,(ci_g(4,2)-ci_g(4,1))/2]; % the pixel distance between the reported and 95% confidence bound for G1 and G2

%% Reporting the centroids
centroids = [fitvariables(2),fitvariables(4)];

%% Calculating the row SNR for both gates
    signal = fitvariables(3);
    SNR = signal/noise;

%% Fitting variables to send out
    fitvariables_out = [fitvariables(1:2),widths(1),fitvariables(3:4),widths(2)];

end

