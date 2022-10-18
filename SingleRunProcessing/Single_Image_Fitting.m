function [centroids,snr,centroid_error,R2,fitvariables,residual,nearwall_bounds,...
    near_wall_extrap,signal] = Single_Image_Fitting(bkg_subtracted_imageData_ROI,...
    cutoff_height_pixels,emissionlocatingdata,gate1_location_bounds,gate2_location_bounds,...
                                    fitting_limits,Delays,Gates,run,resolution,rows,cfd_turb_prof,...
                                    synth_switch)
%This function fits a single image using a bunch of pre-processing
%information. The program works by first fitting data above the wall where
%the gates are far apart. The program then extrapolates the results down
%closer to the wall for the first gate. Finally, the displaced emissions
%are fit using the extrapolated results. 

%%%%%%%%% Initializing variables %%%%%%%%%%%%%%%%
        centroids = zeros(rows,2);
        centroid_error = zeros(rows,2);
        R2 = zeros(rows,1);
        snr = zeros(rows,1);
        signal= zeros(rows,1);
        fitvariables = zeros(rows,6);
        imageData_ROI_backup = bkg_subtracted_imageData_ROI;
        residual = zeros(size(bkg_subtracted_imageData_ROI));

%%%%%%%%%%%%%%%%%% Curve Fitting above the overlap %%%%%%%%%%%%%%%%%%%%%%
            %setup the bounds for near wall and not near wall indexes
            wall_yloc = ceil(emissionlocatingdata(2));
            max_len = size(bkg_subtracted_imageData_ROI,1);
            nearwall_bounds = [wall_yloc-round(cutoff_height_pixels),max_len];
            interp_bounds = [nearwall_bounds(1)-120,nearwall_bounds(1)-1];
            not_wall_vec = [1:nearwall_bounds(1),nearwall_bounds(2):rows];
        
            %fit the rows of data that aren't near the wall (minimal gate overlap)
            for i = not_wall_vec
                fitting_limits(:,2) = gate1_location_bounds(i,:);
                fitting_limits(:,5) = gate2_location_bounds(i,:);
                [centroids(i,:),snr(i),centroid_error(i,:),R2(i),fitvariables(i,:),signal(i,:),residual(i,:)] = Curvefitting_nonlinlsqr_gauss(bkg_subtracted_imageData_ROI(i,:),...
                    fitting_limits);
            end

%%%%%%%%%%%%%% Extrapolating results up to the wall   %%%%%%%%%%%%%%%%%%%%%%%%
%% Fitting gate 1 as a line (substantially above the wall) where there is good SNR
        [near_wall_extrap,g1_centroid_error_extrap,centroids_g1] = NearWallExtrap(interp_bounds,nearwall_bounds,fitvariables,Delays,...
            Gates,resolution,wall_yloc,run,snr,rows,centroid_error,cfd_turb_prof,synth_switch); %provides the extrapolated g1 variables for the set of points contained by the nearwall_bounds
        centroids(:,1) = centroids_g1;

%%%%%%%%%%%%%%%%%% Fitting Below the Switchover   %%%%%%%%%    
extrap_fitting_limits = fitting_limits;
min_amp_g1 = min(near_wall_extrap(interp_bounds(1):interp_bounds(2),1));
     for i = nearwall_bounds(1):nearwall_bounds(2)
            extrap_fitting_limits(:,2:3) = repmat(near_wall_extrap(i,2:3),3,1);
            extrap_fitting_limits(:,5) = gate2_location_bounds(i,:);
            extrap_fitting_limits(:,1) = near_wall_extrap(i,1);

            [centroids(i,:),snr(i),centroid_error(i,:),R2(i),fitvariables(i,:),signal(i,:),residual(i,:)] = Curvefitting_nonlinlsqr_gauss(bkg_subtracted_imageData_ROI(i,:),...
                extrap_fitting_limits);
     end

  centroid_error(:,1) = g1_centroid_error_extrap;

%% Some Clean-up
R2(R2<0) = 0;

end