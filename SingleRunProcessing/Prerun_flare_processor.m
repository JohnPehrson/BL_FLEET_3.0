function [flare_isolated_filtered_out,imageData_ROI_flaresubtracted] = Prerun_flare_processor(prerunData_mean,fitting_limits,emissionlocatingdata,gate1_location_bounds,run,imageData_ROI)
%This function takes the time-averaged data prior to the run alongside some
%fitting bounds and uses it to isolate the reflected light coming out of
%the beam port. 

%% Fitting data
    offset_from_surf = 25;
    rows_process = (emissionlocatingdata(2)-offset_from_surf);
    rowcounter = 1:size(prerunData_mean,2);


    options = optimset(@lsqcurvefit);
    options.Display = 'off';
    options.TolFun = 1e-5;
    options.MaxFunEvals = 1e3;
    options.MaxIter = 200;
    options.FinDiffType = 'central';
    fit_gauss2 = @(p) (p(1)*exp(-(rowcounter-p(2)).^2 ./ (2*p(3)^2)));

    LB = [fitting_limits(1,1:3)];
    x0 = [fitting_limits(2,1:3)];
    UB = [fitting_limits(3,1:3)];
    LB(2) = gate1_location_bounds(1,1);
    x0(2) = gate1_location_bounds(1,2);
    UB(2) = gate1_location_bounds(1,3);

    fitvariables_rows = zeros(size(prerunData_mean,1),3);
    fit_rows = zeros(size(prerunData_mean,1),size(prerunData_mean,2));
for i = 1:rows_process
    onerow_data = prerunData_mean(i,:);
    % Fitting Function
    err_fit_gauss2 = @(v) fit_gauss2(v)-onerow_data;
    % Fitting Process
    [fitvariables_rows(i,:),~,residual,~,~,~,jacobian] = lsqnonlin(err_fit_gauss2,x0,LB,UB,options);
    fit_rows(i,:) = fit_gauss2(fitvariables_rows(i,:));

%         figure;
%         scatter(rowcounter,onerow_data,'k','Linewidth',2);
%         hold on;
%         plot(rowcounter,fit_gauss2(fitvariables_rows(i,:)),':k','Linewidth',2);
%         grid on;
%         set(gca,'FontSize', 20);
%         set(gca,'fontname','times')  % Set it to times
%         disp('wait');

end

% 
% figure;
% title('fit');
% image(fit_rows)
% colorbar;
% colormap(turbo(1500));
% axis equal;
% set(gca, 'YDir','reverse')
% 
% figure;
% title('subt.');
% image(prerunData_mean-fit_rows)
% colorbar;
% colormap(turbo(1500));
% axis equal;
% set(gca, 'YDir','reverse')

%% Extrapolating data

extrap_rows = (rows_process+1):size(prerunData_mean,1);

extrap_height = 80;
prefit_fitrows = (rows_process-extrap_height):rows_process;
p_loc = polyfit(prefit_fitrows,fitvariables_rows(prefit_fitrows,2),1);
p_width = polyfit(prefit_fitrows,fitvariables_rows(prefit_fitrows,3),1);
p_amp = polyfit(prefit_fitrows,fitvariables_rows(prefit_fitrows,1),1);
y_loc = polyval(p_loc,extrap_rows);
y_width = polyval(p_width,extrap_rows);
y_amp = polyval(p_amp,extrap_rows);

test = (abs(length(extrap_rows)-(extrap_rows-min(extrap_rows))))./length(extrap_rows);

        for i = extrap_rows
        if i<(emissionlocatingdata(2)+3)
            scale = 1;
        else
            scale = test(i-rows_process);
        end

        fitvariables_rows(i,:) = [scale.*y_amp(i-rows_process),y_loc(i-rows_process),y_width(i-rows_process)];
        fit_rows(i,:) = fit_gauss2(fitvariables_rows(i,:));
        end
% 
% figure;
% title('ext. subt.');
% image(prerunData_mean-fit_rows)
% colorbar;
% colormap(turbo(3000));
% axis equal;
% set(gca, 'YDir','reverse')

%% Do some filtering to just get the flare
flare_isolated = prerunData_mean-fit_rows;
passdata = (emissionlocatingdata(2)-offset_from_surf):(emissionlocatingdata(2)+offset_from_surf);
flare_isolated_filtered = zeros(size(flare_isolated));
flare_isolated_filtered(passdata,:) = flare_isolated(passdata,:);
flare_isolated_filtered(flare_isolated<100) = 0;

%% Gaussian Blur
flare_isolated_filtered = imgaussfilt(flare_isolated_filtered,2);

% figure;
% title('Flare Surf');
% image(flare_isolated_filtered)
% colorbar;
% colormap(turbo(max(flare_isolated_filtered(:))));
% axis equal;
% set(gca, 'YDir','reverse')

%% Move the flare to where it is during the actual run
if (run==1)
    left_move = 0;
    move_up = 0;
elseif(run==2)
    left_move = 12;
    move_up = 0;
else
    left_move = 0;
    move_up = 0;
end

flare_isolated_filtered_moved = imtranslate(flare_isolated_filtered,[-left_move, -emissionlocatingdata(3)-move_up],'FillValues',0);

% figure;
% title('Flare Surf');
% image(flare_isolated_filtered_moved)
% colorbar;
% colormap(turbo(max(flare_isolated_filtered_moved(:))));
% axis equal;
% set(gca, 'YDir','reverse')

%% Subtraction

if (run==1)
    scale_flare = 0.8;
elseif (run==2)
    scale_flare = 0.4;
else
    scale_flare = 0.4;
end
flare_isolated_filtered_out = scale_flare.*flare_isolated_filtered_moved;
imageData_ROI_flaresubtracted = imageData_ROI-flare_isolated_filtered_out;
imageData_ROI_flaresubtracted(imageData_ROI_flaresubtracted<0) = 0;

figure;
title('Flare Surf');
image(imageData_ROI_flaresubtracted)
colorbar;
colormap(turbo(max(imageData_ROI_flaresubtracted(:))));
axis equal;
set(gca, 'YDir','reverse')

end