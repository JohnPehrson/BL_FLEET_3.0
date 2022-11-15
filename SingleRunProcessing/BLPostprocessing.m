function [mean_velocity,mean_velocity_r,mean_velocity_s,rms_velocity,rms_velocity_r,...
    rms_velocity_s,mean_SNR,row_mm,...
    sufficient_counter,filt_centroids] = BLPostprocessing(red_centroids,...
    red_velocity,red_velocity_s,red_velocity_r,red_g2SNR,...
    gate1_location_bounds,gate2_location_bounds,invar_image_compare,synth_switch,...
    emissionlocatingdata,pixel_um_resolution,...
    y_mm)
%This function postprocesses the fit data. This is broadly devided into
%filtering and calculating time-averaged velocimetry (mean and rms)

    row_mm = y_mm;
%% Filtering the raw fit data
[filt_centroids,filt_velocity,filt_velocity_s,filt_velocity_r,...
    filt_g2SNR,images_percentage_passed_filtering,filt_binary_master] = Filtering(red_centroids,...
    red_velocity,red_velocity_s, red_velocity_r,red_g2SNR,...
    gate1_location_bounds,gate2_location_bounds,synth_switch,row_mm);

%% Calculating the time-averaged data (mean and rms velocity calculations
 [mean_velocity,rms_velocity,mean_velocity_r,mean_velocity_s,rms_velocity_r,rms_velocity_s,...
    mean_SNR,sufficient_counter] = Mean_RMS_Velo_Calculator(filt_velocity,filt_velocity_r,filt_velocity_s,...
    filt_g2SNR,invar_image_compare,filt_binary_master,emissionlocatingdata);

end