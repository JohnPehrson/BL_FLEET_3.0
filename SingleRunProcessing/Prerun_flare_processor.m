function [flare_mask_reg,imageData_ROI_flaresubtracted] = Prerun_flare_processor(prerunData_mean,...
    fitting_limits,emissionlocatingdata,gate1_location_bounds,imageData_mean,flare_scale)
%This function takes the time-averaged data prior to the run alongside some
%fitting bounds and uses it to isolate the reflected light coming out of
%the beam port. 

%% Fitting prerun data
    offset_from_surf = 20;
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

%         figure(3);
%         scatter(rowcounter,onerow_data,'k','Linewidth',2);
%         hold on;
%         plot(rowcounter,fit_gauss2(fitvariables_rows(i,:)),':k','Linewidth',2);
%         grid on;
%         set(gca,'FontSize', 20);
%         set(gca,'fontname','times')  % Set it to times
%         hold off;
end

% close all;
% figure;
% image(fit_rows)
% colorbar;
% colormap(turbo(1500));
% axis equal;
% set(gca, 'YDir','reverse')
% title('fit');
% 
% 
% figure;
% image(prerunData_mean-fit_rows)
% colorbar;
% colormap(turbo(4096));
% axis equal;
% set(gca, 'YDir','reverse')
% title('subt.');
% 

%% Do some filtering to just get the flare
flare_isolated = prerunData_mean-fit_rows;
passdata = (emissionlocatingdata(2)-offset_from_surf):(size(flare_isolated,1));
flare_isolated_filtered = zeros(size(flare_isolated));
flare_isolated_filtered(passdata,:) = flare_isolated(passdata,:);
flare_isolated_filtered(flare_isolated<150) = 0;

%% Locate the Ellipse template over the real image, looking to extract the flare near the wall
half_height = 8; %pix
half_width = 40; %pix

f1 = figure;
image(imageData_mean);
set(gca,'YDir','reverse');
colorbar;
colormap(jet(round(max(imageData_mean(:)))));
hold on;
grid on;
title('Image Mean with ROI')

h = drawellipse('Center',[emissionlocatingdata(1),emissionlocatingdata(2)],'SemiAxes',[half_width,half_height],'StripeColor','r');

mask = createMask(h);
mean_mask = imageData_mean.*mask;

figure;
image(mean_mask)
colorbar;
colormap(turbo(max(mean_mask(:))));
axis equal;
set(gca, 'YDir','reverse')
title('Image Mean Mask')


%% Locate the Ellipse template over the flare image

[~,height_max_loc] = max(sum(flare_isolated_filtered,2)); %for the vertical center location

[amp,span_max_loc] = max(sum(flare_isolated_filtered,1));
x = rowcounter;
x0 = [amp,span_max_loc,half_width];
LB = [200,min(x),half_width/3];
UB = [4096,max(x),2*half_width];
y = sum(flare_isolated_filtered,1);
[fitvariables] = SingleGaussFit(x0,LB,UB,x,y);
span_max_loc = fitvariables(2);

f2 = figure;
image(prerunData_mean);
set(gca,'YDir','reverse');
colorbar;
colormap(jet(round(max(imageData_mean(:)))));
title('Image Prerun Mask')

h = drawellipse('Center',[span_max_loc,height_max_loc],'SemiAxes',[half_width,half_height],'StripeColor','r');


mask = createMask(h);
flare_mask = prerunData_mean.*mask;

figure;
image(flare_mask)
colorbar;
colormap(turbo(max(flare_mask(:))));
axis equal;
set(gca, 'YDir','reverse');
title('Prerun Masked');


%% Image Registration
[optimizer, metric] = imregconfig('monomodal');
flare_mask_reg = imregister(flare_mask,mean_mask,'translation',optimizer,metric);
tform = imregtform(flare_mask,mean_mask,'translation',optimizer,metric);
[x_trans,y_trans] = transformPointsForward(tform,height_max_loc,span_max_loc);
    %manual additional translate to offset for g1-g2 biasing the flare
    %emissions up
    down_force = flare_scale(1); %pixels
    flare_mask_reg = imtranslate(flare_mask_reg,[flare_scale(3), down_force]);
x_rel = x_trans-height_max_loc+down_force; %pos is down
y_rel = y_trans-span_max_loc; %pos is right

figure;
imshowpair(flare_mask_reg,mean_mask);
title('Registration')

%% Scaling the flare to match the time-averaged mean flare emissions
figure;
subplot(1,3,1);
image(flare_mask_reg)
colorbar;
colormap(turbo(max(mean_mask(:))));
grid on;
title('Prerun Mean')
subplot(1,3,2);
image(mean_mask)
colorbar;
colormap(turbo(max(mean_mask(:))));
grid on;
title('Run Mean')


pix_diff = (flare_mask_reg-mean_mask)./flare_mask_reg;
non_nan_pix_diff = pix_diff(~isnan(pix_diff));
scale_down = prctile(non_nan_pix_diff,25);
flare_mask_reg_scale = scale_down.*flare_mask_reg;

subplot(1,3,3);
image(flare_mask_reg_scale)
colorbar;
colormap(turbo(max(mean_mask(:))));
grid on;
title('Prerun Mean, Scaled')

sig = 1;
mean_mask_g = imgaussfilt(mean_mask,sig);
flare_mask_g = imgaussfilt(flare_mask_reg_scale,sig);

figure;
mean_diff = mean_mask_g-flare_mask_g;
image(mean_diff);
colorbar;
colormap(turbo(max(mean_mask_g(:))));
grid on;
title('Flare Subtracted')



%% Subtraction
imageData_ROI_flaresubtracted = imageData_mean-flare_mask_g;
imageData_ROI_flaresubtracted(imageData_ROI_flaresubtracted<0) = 0;

figure;
title('Time-average emissions with the light from the hole subtracted');
image(imageData_ROI_flaresubtracted)
colorbar;
colormap(turbo(max(imageData_ROI_flaresubtracted(:))));
axis equal;
set(gca, 'YDir','reverse')

close([f1 f2])

end