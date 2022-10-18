function [centroids,SNR,centroid_error,R2,fitvariables,total_sig,residual] = Curvefitting_nonlinlsqr_gauss(onerow_data,...
    gaussfit_limits)
%This function fits gaussian curves to a 2-D set of intensity values over a
%pixel range. This function utilizes the nonlinear least squares fitting
%method.

%Just doing 2 gaussian curves right now for G1 and G2.

%% Trim to only use the central part of the image for the curve fitting
rowcounter = 1:length(onerow_data);
blank_row = zeros(1,length(onerow_data));

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
LB = [LB(1:2);LB(3:4)];
x0 = [gaussfit_limits(2,notwidths)];
x0 = [x0(1:2);x0(3:4)];
UB = [gaussfit_limits(3,notwidths)];
UB = [UB(1:2);UB(3:4)];
widths = gaussfit_limits(2,[3,6])';

%% Fitting G1 and G2
        i=1; %gate 1
        %get data for the fitting
        g1_cols = round((x0(i,2)-2*widths(i)):(x0(i,2)+2*widths(i)));
        g1_onerow = onerow_data(g1_cols);
        %fitting functions
        fit_gauss_g1 = @(p) (p(1)*exp(-(g1_cols-p(2)).^2 ./ (2*widths(i)^2)));
        err_fit_gauss_g1 = @(v) fit_gauss_g1(v)-g1_onerow;
        %fit
        [fitvariables_g1,~,residual_g1,~,~,~,jacobian_g1] = lsqnonlin(err_fit_gauss_g1,x0(i,:),LB(i,:),UB(i,:),options);
        %subtract the fit from the rest of data processing
        placed_fit = blank_row;
        placed_fit(g1_cols) = fit_gauss_g1(fitvariables_g1);
        onerow_data = onerow_data-placed_fit;

        i=2; %gate 2
       %get data for the fitting
        g2_cols = round(LB(i,2):UB(i,2));
        g2_onerow = onerow_data(g2_cols);
        %fitting functions
        fit_gauss_g2 = @(p) (p(1)*exp(-(g2_cols-p(2)).^2 ./ (2*widths(i)^2)));
        err_fit_gauss_g2 = @(v) fit_gauss_g2(v)-g2_onerow;
        %fit
        [fitvariables_g2,~,residual_g2,~,~,~,jacobian_g2] = lsqnonlin(err_fit_gauss_g2,x0(i,:),LB(i,:),UB(i,:),options);


%% Confidence in fit
%centroids
ci_g1 = nlparci(fitvariables_g1,residual_g1,'jacobian',jacobian_g1);   %  95% confidence intervals for the fit coefficients
ci_g2 = nlparci(fitvariables_g2,residual_g2,'jacobian',jacobian_g2);   %  95% confidence intervals for the fit coefficients
centroid_error = [(ci_g1(2,2)-ci_g1(2,1))/2,(ci_g2(2,2)-ci_g2(2,1))/2]; % the pixel distance between the reported and 95% confidence bound for G1 and G2

%% Reporting the centroids
centroids = [fitvariables_g1(2),fitvariables_g2(2)];

%% Calculating the row SNR for both gates
    z_90 = 1.645; %capture the central 90% of the 2nd gate
    signal_vec_g2 = fit_gauss_g2(fitvariables_g2);
    fit_residual_vec_g2 = abs(err_fit_gauss_g2(fitvariables_g2));
    signal = trapz(signal_vec_g2);
    noise = trapz(fit_residual_vec_g2);
    SNR = signal/noise;

    total_sig = sum(signal_vec_g2);
    R2 = 0;

    %% Residual
    r_g1 = blank_row;
    r_g2 = blank_row;
    r_g1(g1_cols) = residual_g1;
    r_g2(g2_cols) = residual_g2;
    residual = r_g1+r_g2;

%% Put the fit variables back in the right order
fitvariables = [fitvariables_g1,widths(1),fitvariables_g2,widths(2)];

end

