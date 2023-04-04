function [centroids,snr,centroid_error,fitvariables,residual,nearwall_bounds,...
    near_wall_extrap,difference_imfit] = Single_Image_Fitting(bkg_subtracted_imageData_ROI,...
    cutoff_height_pixels,emissionlocatingdata,gate1_location_bounds,gate2_location_bounds,...
                                    fitting_limits,Delays,Gates,run,resolution,rows,cfd_turb_prof,...
                                    synth_switch,near_wall_g1_scale,noise)
%This function fits a single image using a bunch of pre-processing
%information. The program works by first fitting data above the wall where
%the gates are far apart. The program then extrapolates the results down
%closer to the wall for the first gate. Finally, the displaced emissions
%are fit using the extrapolated results. 

%%%%%%%%% Initializing variables %%%%%%%%%%%%%%%%
        centroids = zeros(rows,2);
        centroid_error = zeros(rows,2);
        snr = zeros(rows,1);
        fitvariables = zeros(rows,6);
        residual = zeros(size(bkg_subtracted_imageData_ROI));

%%%%%%%%%%%%%%%%%% Curve Fitting above the overlap %%%%%%%%%%%%%%%%%%%%%%
            %setup the bounds for near wall and not near wall indexes
            wall_yloc = ceil(emissionlocatingdata(2));
            max_len = size(bkg_subtracted_imageData_ROI,1);
            nearwall_bounds = [round(cutoff_height_pixels),max_len];
            interp_bounds = [1,nearwall_bounds(1)-1];
            not_wall_vec = [1:nearwall_bounds(1)];
        
            %fit the rows of data that aren't near the wall (minimal gate overlap)
            for i = not_wall_vec
                fitting_limits(:,2) = gate1_location_bounds(i,:);
                fitting_limits(:,5) = gate2_location_bounds(i,:);
                [centroids(i,:),snr(i),centroid_error(i,:),fitvariables(i,:),residual(i,:)] = Curvefitting_nonlinlsqr_gauss(bkg_subtracted_imageData_ROI(i,:),...
                    fitting_limits,i,noise);
            end

%% Centroid Gate 1 as a line
high_snr_binary = snr>(prctile(snr,50));
c1_high_snr_rows = not_wall_vec(high_snr_binary);
c1 = fitvariables(c1_high_snr_rows,2);
p = polyfit(c1_high_snr_rows,c1,1);
allrows = 1:rows;
c1_line = transpose(polyval(p,allrows));
fitvariables(:,2) = c1_line;

%% Extrapolating g1 location and amplitude up to the wall
        [near_wall_extrap,g1_centroid_error_extrap,centroids_g1] = NearWallExtrap(interp_bounds,nearwall_bounds,fitvariables,Delays,...
            Gates,resolution,wall_yloc,run,snr,rows,centroid_error,cfd_turb_prof,synth_switch,near_wall_g1_scale); %provides the extrapolated g1 variables for the set of points contained by the nearwall_bounds
        centroids(:,1) = centroids_g1;

%% Plotting the image with the fit and extrapolated fit to g1 subtracted. This gives only the data that is being fit for the second gate
cols_list = 1:size(bkg_subtracted_imageData_ROI,2);
rows_list = 1:size(bkg_subtracted_imageData_ROI,1);
g1_fit = zeros(size(bkg_subtracted_imageData_ROI));
fit_gauss = @(p) (p(1)*exp(-(cols_list-p(2)).^2 ./ (2*p(3)^2)));
for i = 1:length(rows_list)
    g1_row_fit = fit_gauss(near_wall_extrap(i,:));
    g1_fit(i,:) = g1_row_fit;
end

difference_imfit = bkg_subtracted_imageData_ROI-g1_fit;
% 
% figure;
% subplot(1,2,1);
% image(g1_fit);
% colorbar;
% colormap(jet(round(max(g1_fit(:)))));
% subplot(1,2,2);
% image(difference_imfit);
% colorbar;
% colormap(jet(round(max(g1_fit(:)))));

%%%%%%%%%%%%%%%%%% Fitting Below the Switchover   %%%%%%%%%    
extrap_fitting_limits = fitting_limits;
     for i = nearwall_bounds(1):nearwall_bounds(2)
            extrap_fitting_limits(:,2:3) = repmat(near_wall_extrap(i,2:3),3,1);
            extrap_fitting_limits(:,5) = gate2_location_bounds(i,:);
            extrap_fitting_limits(:,1) = near_wall_extrap(i,1);

            [centroids(i,:),snr(i),centroid_error(i,:),fitvariables(i,:),residual(i,:)] = Curvefitting_nonlinlsqr_gauss(bkg_subtracted_imageData_ROI(i,:),...
                extrap_fitting_limits,i,noise);
     end

  centroid_error(:,1) = g1_centroid_error_extrap;
end