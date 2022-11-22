function [red_centroids,red_velocity,red_velocity_s,red_velocity_r,red_g2SNR,...
    uncertainty_wall_loc_pix,numimages] = BLImageProcessor(run,filepath,...
    ROI,imageData_mean,dust_filter_bounds_bottom,dust_filter_bounds_top,gate1_location_bounds,...
    gate2_location_bounds,background_totalfit,emissionlocatingdata,cutoff_height_pixels,pixel_um_resolution,gates,...
    delay,fitting_limits,synthfilepath,synth_switch,constant_background,rotation_angles,...
    total_sig_thresh,cfd_turb_prof,imageprocess_numbers,near_wall_g1_scale)
%This function performs the time-accurate FLEET measurements by gaussian
%fitting each row of each image to find the location of the emissions at
%specific times. These locations are used to find instantaneous velocities.
cols = ROI(4)-ROI(3)+1;
rows = ROI(2)-ROI(1)+1;

%% Initializing reduction variables for image processing
%vectors
red_centroids = []; %centroid locations as f(y) for each image
red_velocity = [];  %velocity as f(y) for each image
red_velocity_s = [];  %systematic velocity as f(y) for each image
red_velocity_r = [];  %random velocity as f(y) for each image
red_R2 = [];        %R2 as f(y) for each image
red_g2SNR = [];  %SNR as f(y) for each image
red_signal =[]; %signal as f(y) for each image
red_move_ROI_ud = []; %correction in pixels

    %% Parallel processing
    WaitMessage = parfor_wait(length(imageprocess_numbers), 'Waitbar', true);
    parfor imloop = 1:length(imageprocess_numbers)
        
        %Initializing temporary variables
        %just used locally
            temp_centroids = zeros(2,rows);
            temp_velocity = zeros(2,rows);
            temp_velocity_s = zeros(2,rows);
            temp_velocity_r = zeros(2,rows);
            temp_R2 = zeros(2,rows);
            temp_snr = zeros(2,rows);
            temp_signal = zeros(1,rows);
            temp_move_ROI_ud = zeros(1,rows);

        %Load image
            imagenumber = imageprocess_numbers(imloop);
            [imageData] = ImageLoader(filepath(1),filepath(2),synthfilepath,synth_switch,imagenumber);
            if synth_switch
            imageData_ROI = imageData;
            else 
            imageData = imrotate(imageData,rotation_angles,'bilinear','crop');
            imageData_ROI = imageData(ROI(1):ROI(2),ROI(3):ROI(4)); %get ROI
            end
    %filtering
    total_roi_sig = sum(imageData_ROI(:));
        if and(total_roi_sig>total_sig_thresh(1),total_roi_sig<total_sig_thresh(2))
             %Image Registration (moving the image
            [imageData_ROI,temp_move_ROI_ud] = ImageRegistration(imageData_ROI,imageData_mean,emissionlocatingdata);
    
            [noise] = NoiseCalculator(imageData_ROI);

            %Filter image and subtract background
            imageData_ROI(imageData_ROI<dust_filter_bounds_bottom) = imageData_mean(imageData_ROI<dust_filter_bounds_bottom);
            imageData_ROI(imageData_ROI>dust_filter_bounds_top) = imageData_mean(imageData_ROI>dust_filter_bounds_top);
            imageData_ROI = imageData_ROI-background_totalfit;
            imageData_ROI(imageData_ROI<0) = 0;
    
            %Process image
            [temp_centroids,temp_velocity,temp_velocity_r,temp_velocity_s,...
                temp_snr] = BLImageCurveFitter_VelocityCalculator(imageData_ROI,cols,...
                rows,gate1_location_bounds,gate2_location_bounds,emissionlocatingdata,...
                fitting_limits,cutoff_height_pixels,pixel_um_resolution,gates,delay,...
                constant_background,run,cfd_turb_prof,synth_switch,near_wall_g1_scale,noise);
    
            %% Save Image and centroids for video
%             SaveImage(run,imageprocess_numbers,imloop,imageData_ROI,temp_centroids)
    
            %% Sort temporary variables into the reduction variables
            red_centroids = [red_centroids,temp_centroids]; %centroid locations as f(y) for each image
            red_velocity = [red_velocity,temp_velocity];    %velocity as f(y) for each image
            red_velocity_s = [red_velocity_s,temp_velocity_s];    %velocity as f(y) for each image
            red_velocity_r = [red_velocity_r,temp_velocity_r];    %velocity as f(y) for each image
            red_g2SNR = [red_g2SNR,temp_snr];      %SNR as f(y) for each image
            red_move_ROI_ud = [red_move_ROI_ud,temp_move_ROI_ud];
    
        end
       WaitMessage.Send;
    end
    WaitMessage.Destroy
    uncertainty_wall_loc_pix = std(red_move_ROI_ud);
    numimages = size(red_velocity,2);

end