function [g1_location_col,wallfit_location_col,move_flare_fit,zero_height_ref_unc] = MeanLocational_Fitting(imageData_mean,...
    drow,dcol,run,prerunData_mean,top_offset)
%This function is meant to get time-averaged results about where the FLEET
%line lies. This function returns data that allows instantaneous FLEET
%measurements to be spatially located consistently in the x and y axis.

%% Locate emissions in dx (finding the location of the first gate) 
colaverage = mean(imageData_mean,1);
% figure;
% plot(1:dcol,colaverage);

[~,maxind] = max(colaverage);
g1_location_col = maxind;


%% Location the emissions in dy (finding the location of the bright spot near the surface) for the real data

    %max
        x = 1:drow;
        sum_half_width = 12;
        rowaverage = mean(imageData_mean(:,round(g1_location_col-sum_half_width):round(g1_location_col+sum_half_width)),2);
        TF_row = islocalmax(rowaverage,'MinProminence',10);
        possible_wall_locations = x(TF_row);
        [~,closest_wall_loc_ind] = min(abs(top_offset-possible_wall_locations));
        wallfit_location_col = possible_wall_locations(closest_wall_loc_ind);
        amp = rowaverage(wallfit_location_col);

        figure;
        plot(x,rowaverage);
        hold on;
        scatter(x(TF_row),rowaverage(TF_row))

    %gaussian fit
        x0 = [amp,wallfit_location_col,2];
        LB = [0.5.*amp,wallfit_location_col-5,1];
        UB = [2.*amp,wallfit_location_col+5,3];
     
            options = optimset(@lsqcurvefit);
            options.Display = 'off';
            options.TolFun = 1e-5;
            options.MaxFunEvals = 1e3;
            options.MaxIter = 200;
            options.FinDiffType = 'central';
            
            fit_gauss_g1 = @(p) (p(1)*exp(-(x-p(2)).^2./(2*p(3)^2)));
            err_fit_gauss_g1 = @(v) fit_gauss_g1(v)-rowaverage';
    
            [fitvariables,~,residual,~,~,~,jacobian] = lsqnonlin(err_fit_gauss_g1,x0,LB,UB,options);
            ci_g1 = nlparci(fitvariables,residual,'jacobian',jacobian);   %  95% confidence intervals for the fit coefficients
            zero_height_ref_unc = [(ci_g1(2,2)-ci_g1(2,1))/2]; 

%% Location the emissions in dy (finding the location of the bright spot near the surface) for the prerun data
   
    %max
        sum_half_width = 12;
        rowaverage = mean(prerunData_mean(:,round(g1_location_col-sum_half_width):round(g1_location_col+sum_half_width)),2);
        [amp2,wallfit_location_col2] = max(rowaverage);

%     %gaussian fit
%         x0 = [amp2,wallfit_location_col2,2];
%         LB = [0.5.*amp2,wallfit_location_col2-5,1];
%         UB = [2.*amp2,wallfit_location_col2+5,3];
%         [fitvariables2] = SingleGaussFit(x0,LB,UB,x,rowaverage');
%         
%         figure;
%         plot(x,rowaverage');
%         hold on;
%         plot(x,fit_gauss(fitvariables2));

%% Difference in height from background to real run
move_flare_fit = wallfit_location_col2-wallfit_location_col;%move up by this amount
end