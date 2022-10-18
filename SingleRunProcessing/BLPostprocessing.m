function [mean_velocity,mean_velocity_r,mean_velocity_s,rms_velocity,rms_velocity_r,...
    rms_velocity_s,mean_SNR,mean_R2,mean_signal,row_mm,row_mm_uncertainty,...
    sufficient_counter,filt_centroids] = BLPostprocessing(red_centroids,...
    red_velocity,red_velocity_s,red_velocity_r,red_R2,red_g2SNR,red_signal,...
    gate1_location_bounds,gate2_location_bounds,invar_image_compare,synth_switch,...
    inclination_angle,inclination_angle_unc,emissionlocatingdata,pixel_um_resolution,zero_height_ref_unc)
%This function postprocesses the fit data. This is broadly devided into
%filtering and calculating time-averaged velocimetry (mean and rms)

if (~synth_switch) %only save if real data
temp_savename = 'Matfiles_preprocessing/timeaccurateunfiltereddata.mat';
save(temp_savename);
end

row_pixels = (1:size(red_centroids,1))';
row_mm = (emissionlocatingdata(2)-row_pixels+1)*pixel_um_resolution(1)./(1000);

%% Unfiltered Plotting 
% g1_centroids = red_centroids(:,1:2:end);
% g2_centroids = red_centroids(:,2:2:end);
% unfiltered_velo_mean = mean(red_velocity,2);
% unfiltered_red_velocity_r = mean(red_velocity_r,2);
% 
% figure;
% plot(unfiltered_velo_mean,row_mm);
% title('Unfiltered Mean velocity')
% 
% figure;
% plot(unfiltered_red_velocity_r,row_mm);
% title('Unfiltered Mean velocity uncertainty')
% 
% figure;
% plot(std(red_velocity'),row_mm);
% title('Unfiltered RMS velocity')


%% Filtering the raw fit data
[filt_centroids,filt_velocity,filt_velocity_s,filt_velocity_r,filt_R2,...
    filt_g2SNR,filt_signal,images_percentage_passed_filtering,filt_binary_master] = Filtering(red_centroids,...
    red_velocity,red_velocity_s, red_velocity_r,red_R2,red_g2SNR,red_signal,...
    gate1_location_bounds,gate2_location_bounds,synth_switch,row_mm);

    %% Save data for time-accurate velocimetry (for comparison against the DAQ)
    rows = size(filt_velocity,1);
    above_surf_rows = (round(emissionlocatingdata(2)):-1:1).*(pixel_um_resolution(1)./1000);
    frstream_rows = (above_surf_rows<12)&(above_surf_rows>8);
    time_accurate_freestream_velo = zeros(1,size(filt_velocity,2));
    
    for i = 1:length(time_accurate_freestream_velo)
            above_surf = filt_velocity(1:length(above_surf_rows),i);
            freestreamvelocity = above_surf(frstream_rows);
            freestreamvelocity = freestreamvelocity(~isnan(freestreamvelocity));
            if length(freestreamvelocity)<3
            freestreamvelocity = NaN;
            end
            time_accurate_freestream_velo(i) = mean(freestreamvelocity);
    end

    if length(time_accurate_freestream_velo)>500
    savename = "C:\Users\clark\Documents\GitHub\BL_FLEET\SyntheticDataCorrection\TimeAccurateFLEET\timeaccuratedata.mat";
    save(savename,'time_accurate_freestream_velo')
    end

%% Calculating the time-averaged data (mean and rms velocity calculations
 [mean_velocity,rms_velocity,mean_velocity_r,mean_velocity_s,rms_velocity_r,rms_velocity_s,...
    mean_SNR,mean_R2,mean_signal,sufficient_counter] = Mean_RMS_Velo_Calculator(filt_velocity,filt_velocity_r,filt_velocity_s,...
    filt_R2,filt_g2SNR,filt_signal,invar_image_compare,filt_binary_master,emissionlocatingdata);

%% Add in camera pointing uncertainty (inclination, panning, wall location)
[row_mm,row_mm_uncertainty,mean_velocity_s,rms_velocity_s] = CameraPointingProp(mean_velocity,...
 rms_velocity,mean_velocity_s,rms_velocity_s,inclination_angle,inclination_angle_unc,...
 emissionlocatingdata,pixel_um_resolution,row_mm,zero_height_ref_unc);

end