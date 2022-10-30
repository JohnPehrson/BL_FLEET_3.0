function [centroids,snr,centroid_error,R2,fitvariables,residual,nearwall_bounds,...
    near_wall_extrap,signal] = Single_Image_Fitting(bkg_subtracted_imageData_ROI,...
    cutoff_height_pixels,emissionlocatingdata,gate1_location_bounds,gate2_location_bounds,...
                                    fitting_limits,Delays,Gates,run,resolution,rows,cfd_turb_prof,...
                                    synth_switch,near_wall_g1_scale)
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
            nearwall_bounds = [round(cutoff_height_pixels),max_len];
            interp_bounds = [1,nearwall_bounds(1)-1];
            not_wall_vec = [1:nearwall_bounds(1)];
        
            %fit the rows of data that aren't near the wall (minimal gate overlap)
            for i = not_wall_vec
                fitting_limits(:,2) = gate1_location_bounds(i,:);
                fitting_limits(:,5) = gate2_location_bounds(i,:);
                [centroids(i,:),snr(i),centroid_error(i,:),R2(i),fitvariables(i,:),signal(i,:),residual(i,:)] = Curvefitting_nonlinlsqr_gauss(bkg_subtracted_imageData_ROI(i,:),...
                    fitting_limits,i);
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

%%%%%%%%%%%%%%%%%% Fitting Below the Switchover   %%%%%%%%%    
extrap_fitting_limits = fitting_limits;
     for i = nearwall_bounds(1):nearwall_bounds(2)
            extrap_fitting_limits(:,2:3) = repmat(near_wall_extrap(i,2:3),3,1);
            extrap_fitting_limits(:,5) = gate2_location_bounds(i,:);
            extrap_fitting_limits(:,1) = near_wall_extrap(i,1);

            [centroids(i,:),snr(i),centroid_error(i,:),R2(i),fitvariables(i,:),signal(i,:),residual(i,:)] = Curvefitting_nonlinlsqr_gauss(bkg_subtracted_imageData_ROI(i,:),...
                extrap_fitting_limits,i);
     end

  centroid_error(:,1) = g1_centroid_error_extrap;

%% Some Clean-up
R2(R2<0) = 0;

end