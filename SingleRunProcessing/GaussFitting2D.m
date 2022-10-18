function [Data2d_subtracted,fit_values,flare_fit_variables] = GaussFitting2D(Data2d,rows,cols,x0)
%This function should be used to fit the higher surface intensity near the
%surface.

    %% 3D fitting
    %Fitting bounds
    width_wiggle = 1;
    wall_wiggle = 3;
    vertical_wiggle = 2;
    LB = [x0(1),x0(2)-vertical_wiggle,x0(3)-width_wiggle,x0(4)-wall_wiggle,x0(5)-width_wiggle];
    UB = [x0(1),x0(2)+vertical_wiggle,x0(3)+width_wiggle,x0(4)+wall_wiggle,x0(5)+width_wiggle];

    %model for local noise
    [YY,XX] = meshgrid(1:cols,1:rows);
    fit_gauss = @(p) (p(1)*exp(-(XX-p(2)).^2 ./ (2*p(3)^2)).*(exp(-(YY-p(4)).^2 ./ (2*p(5)^2))));
    err_fit_gauss = @(v) fit_gauss(v)-Data2d;

    %Options
    options = optimset(@lsqcurvefit);
    options.Display = 'off';
    options.TolFun = 1e-5;
    options.MaxFunEvals = 1e3;
    options.MaxIter = 200;
    options.FinDiffType = 'central';

    %fit
    [flare_fit_variables] = lsqnonlin(err_fit_gauss,x0,LB,UB,options);
    fit_values = fit_gauss(flare_fit_variables);    

    Data2d_subtracted = Data2d-fit_values;
    Data2d_subtracted(Data2d_subtracted<0) = 0;

    plot_amp = 1500;
    figure;
    subplot(1,3,1);
    image(Data2d)
    colorbar;
    colormap(bone(plot_amp));
    title('Data')

    subplot(1,3,2);
    image(fit_values)
    colorbar;
    colormap(bone(plot_amp));
    title('Fit')

    subplot(1,3,3);
    image(Data2d_subtracted)
    colorbar;
    colormap(bone(plot_amp));
    title('Difference')

    disp('wait');
    
end