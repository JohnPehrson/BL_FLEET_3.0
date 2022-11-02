function [flare_mask_reg_scale,imageData_ROI_flaresubtracted] = Prerun_flare_processor(prerunData_mean,...
    fitting_limits,emissionlocatingdata,gate1_location_bounds,imageData_mean,flare_scale,...
    create_prerun_flare_dataset,single_run)
%This function takes the time-averaged data prior to the run alongside some
%fitting bounds and uses it to isolate the reflected light coming out of
%the beam port. 
        
        %% Settings
        half_height = 10; %pix
        half_width = 50; %pix
        sigma_gfilt = [1.5,3];
        savename = ['Preprocessing_Filestorage/Flare_fit.mat'];

if create_prerun_flare_dataset
        
        %% Fitting prerun data further from the surface
            offset_from_surf = 20;
            rows_process = (emissionlocatingdata(2)-offset_from_surf);
            rowcounter = 1:size(prerunData_mean,2);
        
            options = optimset(@lsqcurvefit);
            options.Display = 'off';
            options.TolFun = 1e-5;
            options.MaxFunEvals = 1e3;
            options.MaxIter = 200;
            options.FinDiffType = 'central';
            fit_gauss1 = @(p) (p(1)*exp(-(rowcounter-p(2)).^2 ./ (2*p(3)^2)));
        
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
            err_fit_gauss2 = @(v) fit_gauss1(v)-onerow_data;
            % Fitting Process
            [fitvariables_rows(i,:),~,residual,~,~,~,jacobian] = lsqnonlin(err_fit_gauss2,x0,LB,UB,options);
            fit_rows(i,:) = fit_gauss1(fitvariables_rows(i,:));
        % 
        %         figure(3);
        %         scatter(rowcounter,onerow_data,'k','Linewidth',2);
        %         hold on;
        %         plot(rowcounter,fit_gauss2(fitvariables_rows(i,:)),':k','Linewidth',2);
        %         grid on;
        %         set(gca,'FontSize', 20);
        %         set(gca,'fontname','times')  % Set it to times
        %         hold off;
        end
        
        %% Extrapolating G1 all the way to the surface 
        use_rows = 50;
        extrap_min = size(prerunData_mean,1);
        extrap_rows = (rows_process+1):extrap_min;
        rows_fit = (rows_process-use_rows):rows_process;
        use_fitvariables_rows = fitvariables_rows(rows_fit,:);
        
            %g1 as a line
                cent_eq = polyfit(rows_fit,use_fitvariables_rows(:,2),1);
                cent_extrap = transpose(polyval(cent_eq,extrap_rows));
            %amp as a line
                amp_eq = polyfit(rows_fit,use_fitvariables_rows(:,1),1);
                amp_extrap = transpose(polyval(amp_eq,extrap_rows));
            %width as a constant
                wid = mean(use_fitvariables_rows(:,3));
                wid_extrap = wid.*ones(length(extrap_rows),1);
        
            %using the extrapolated values to actually fit those rows
            for i = extrap_rows
                alt_i = i-rows_process;
                p = [amp_extrap(alt_i),cent_extrap(alt_i),wid_extrap(alt_i)];
                   fit_rows(i,:) = fit_gauss1(p);
            end
        
            %subtracting the fit from the average to just leave the reflected flare
            %light
            flare_isolated = prerunData_mean-fit_rows;
        
        
%         %% Plotting the prerun flare with and without the fitting
%         close all;
%         figure;
%         subplot(1,2,1);
%         image(prerunData_mean)
%         colorbar;
%         colormap(jet(round(max(prerunData_mean(:)))));
%         axis equal;
%         set(gca, 'YDir','reverse')
%         title('Prerun Image')
%         
%         
%         subplot(1,2,2);
%         image(flare_isolated)
%         colorbar;
%         colormap(jet(round(max(prerunData_mean(:)))));
%         axis equal;
%         set(gca, 'YDir','reverse')
%         title('Prerun Image with G1 subtracted')
        
        %% Isolate the prerun flare using an ellipse
            %basic filtering to get the general area of the flare
                passdata = (emissionlocatingdata(2)-offset_from_surf):(size(flare_isolated,1));
                flare_isolated_filtered = zeros(size(flare_isolated));
                flare_isolated_filtered(passdata,:) = flare_isolated(passdata,:);
                flare_isolated_filtered(flare_isolated<150) = 0;
            %finding the location to place the center of the flare-ellipse
                [amp,span_max_loc] = max(sum(flare_isolated_filtered,1));
                x = rowcounter;
                x0 = [amp,span_max_loc,half_width];
                LB = [200,min(x),half_width/3];
                UB = [4096,max(x),2*half_width];
                y = sum(flare_isolated_filtered,1);
                [fitvariables] = SingleGaussFit(x0,LB,UB,x,y);
                span_max_loc = round(fitvariables(2));
        
                [~,height_max_loc] = max(sum(flare_isolated_filtered,2)); %for the vertical center location
        
        f2 = figure;
        image(flare_isolated);
        set(gca,'YDir','reverse');
        colorbar;
        colormap(jet(round(max(prerunData_mean(:)))));
        hold on;
        grid on;
        title('Prerun Image with G1 subtracted')
        
        h = drawellipse('Center',[span_max_loc,height_max_loc],'SemiAxes',[half_width,half_height],'StripeColor','r');
        mask = createMask(h);
        prerun_mask = flare_isolated.*mask;
        
        %trim mean mask to test sizes
        xbounds = round((span_max_loc-half_width-sigma_gfilt(2)):(span_max_loc+half_width+sigma_gfilt(2)));
        xbounds = xbounds(xbounds>0);
        ybounds = round((height_max_loc-half_height-sigma_gfilt(1)):(height_max_loc+half_height-sigma_gfilt(1)));
        prerun_mask = prerun_mask(ybounds,xbounds);

         %filtering to only get the brightest part of the flare
        prerun_mask_bright = prerun_mask.*(prerun_mask>prctile(prerun_mask(:),40));
        prerun_mask = imgaussfilt(prerun_mask,sigma_gfilt);

%         figure;
%         image(prerun_mask_bright);
%         set(gca,'YDir','reverse');
%         colorbar;
%         colormap(jet(round(max(prerun_mask_bright(:)))));
%         hold on;
%         grid on;
%         title('Prerun Image with G1 subtracted')

        %% Saving
        if single_run==5495
            save(savename,'prerun_mask','single_run');
        end
else %load in data
    load(savename,'prerun_mask');
end

%% Locate the Ellipse template over the real image, looking to extract the flare near the wall
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

%% Comparing the prerun flare to the actual run flare
% figure;
% subplot(1,2,1);
% image(prerun_mask)
% colorbar;
% colormap(turbo(max(prerun_mask(:))));
% axis equal;
% set(gca, 'YDir','reverse')
% title('Image Prerun Mask')
% 
% subplot(1,2,2);
% image(mean_mask)
% colorbar;
% colormap(turbo(max(prerun_mask(:))));
% axis equal;
% set(gca, 'YDir','reverse')
% title('Image Mean Mask')

%% Put the prerun flare fit into a matrix the size of the mean
flare_mat = zeros(size(mean_mask));
prerun_x = size(prerun_mask,2);
prerun_y = size(prerun_mask,1);
mean_x = size(mean_mask,2);
mean_y = size(mean_mask,1);

side_x = round((mean_x-prerun_x)/2);
side_y = round(emissionlocatingdata(2)-(prerun_y/2));
xlocs = side_x:(side_x+prerun_x-1);
ylocs = side_y:(side_y+prerun_y-1);
flare_mat(ylocs,xlocs) = prerun_mask;


%% Image Registration of the mean and template
    [optimizer, metric] = imregconfig('monomodal');
    flare_mask_reg = imregister(flare_mat,mean_mask,'translation',optimizer,metric);
    tform = imregtform(flare_mat,mean_mask,'translation',optimizer,metric);
    [y_trans,x_trans] = transformPointsForward(tform,emissionlocatingdata(2),emissionlocatingdata(1));
    x_rel = x_trans-emissionlocatingdata(1); %pos is down
    y_rel = y_trans-emissionlocatingdata(2); %pos is right

    %manual additional translate to offset for g1-g2 biasing the flare
    %emissions up
    down_force = flare_scale(1); %pixels
    right_force = flare_scale(3); %pixels
    flare_mask_reg = imtranslate(flare_mat,[y_rel+right_force,x_rel+down_force]);

% figure;
% imshowpair(flare_mask_reg,mean_mask);
% title('Registration')

%% Scaling
flare_mask_reg_scale = flare_mask_reg.*flare_scale(2);

%% Subtract template from the time-averaged mean
imageData_ROI_flaresubtracted = imageData_mean-flare_mask_reg_scale;
imageData_ROI_flaresubtracted(imageData_ROI_flaresubtracted<0) = 0;

% figure;
% title('Time-average emissions with the light from the hole subtracted');
% image(imageData_ROI_flaresubtracted)
% colorbar;
% colormap(turbo(max(imageData_ROI_flaresubtracted(:))));
% axis equal;
% set(gca, 'YDir','reverse')

close([f1])

end