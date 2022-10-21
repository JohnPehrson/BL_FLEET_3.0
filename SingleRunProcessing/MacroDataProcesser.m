clear all;close all;clc;

%% Boundary Layer FLEET Analysis Code
    %This code is used for the calculation of mean and fluctuating velocity
    %profiles for Femtosecond Laser Electronic Excitation Tagging (FLEET)
    %taken through a boundary layer. 
    %John C. Pehrson
    %October 20, 2022

%Header

% Navigate to the PreRunConditions Folder and run the 4 scripts in that
% folder
    %Run Conditions - provides IRO conditions, Run labels, FLEET beam
    %location, and the spanwise camera angle reletive the normal

    %Grid Card Reolution Calculator - picks out the run numbers, finds the
    %image resolution and uncertainty, finds the wall location, tabulates
    %the inclination of the camera downward, and identifies the image
    %rotation reletive the freestream

    %ACE on Condition Bounds Processor - provides the reynolds number with
    %uncertainty for the 'steady state' of the run. Provides start and stop
    %images for the run and the pre-run for image processing

    %Filepath Loader loads in an excel file with the filepaths of the FLEET
    %images

%% Input Variables
run = 2;                                                %1: Laminar, 2: 100% Pizza Box, 3: Synthetic Roughness

load('TestConditions/FLEETFilePaths.mat');
run_start_end = [2000,28000;4250,32000;6800,30000]; 
run_start_end_prerun_flare = [300,800;1,1450;1000,4000]; 
runconditions_filepath = 'TestConditions/BLFLEETRunConditions.mat';    %stuff like gates and delays
resolution_filepath = 'TestConditions/RefData.mat';                    %resolution
CFD_turbulent_profile_filepath = 'C:\Users\clark\Documents\GitHub\BL_FLEET_2.0\RawDataProcessing\SynthData\Case056NoTunnel.xlsx';
ROI = [1,1024,300,450;
       250,675,300,440;
       250,675,300,440]; %y limits, then x

fitting_limits = [25,NaN,4,10,NaN,7;...
                  2000,NaN,4.5,1000,NaN,8;...
                  4096,NaN,5,3000,NaN,8]; %LB,x0,UB
total_im_process = 1000; %NaN if process all
total_im_preprocess = 1000; 
loc_data = [140.89,1.00;0,1.00]; %x and uncertainty; y and uncertainty. origin is in the center of the tripping array, x is streamwise
lens_standoff_dist = 415; %mm

if isnan(total_im_process)
total_im_process = run_start_end(run,2)-run_start_end(run,1)-1;
end

    %% Synthetic Data for use if the switch is flipped
    synth_switch = false;  
    synth_numimages = 300;
    if synth_switch
        f = fullfile('ProcessedData',['ProcessedData_Run',num2str(run),'*']);
        fstruct = dir(f);
        synth_real_replicate = fullfile('ProcessedData',fstruct.name);

%         total_im_process(isnan(total_im_process)) = length(run_start_end(run,1):(run_start_end(run,2)-1));
%         synth_real_replicate = ['ProcessedData/ProcessedData_Run',num2str(run),'_Images',num2str(total_im_process),'.mat'];
        total_im_process = synth_numimages;

        [synthfilepath,run_start_end_synth,nondim_velo_error,synth_input_tau_fit,...
            synth_input_velocity_mean] = Synthetic_Data_Gen(runconditions_filepath,...
         resolution_filepath,run,synth_real_replicate,synth_numimages,ROI(run,:),fitting_limits);
        run_start_end(run,:) = run_start_end_synth(run,:);
    else
        synthfilepath = '';
        synth_numimages = total_im_process;
        nondim_velo_error = 0;
        synth_input_tau_fit = 0;
        synth_input_velocity_mean = 0;
        synth_real_replicate = 0;
    end

    %% Image processing (if I need to do it)
    if (~synth_switch)
    data_processing_filepath = ['Matfiles_preprocessing/ProcessedImageData_Run',num2str(run),'Imcount_',num2str(total_im_process),'.mat'];
    else
    data_processing_filepath = ['Matfiles_preprocessing/ProcessedSynthData_Run',num2str(run),'Imcount_',num2str(total_im_process),'.mat'];
    end
    
    if (isfile(data_processing_filepath)) %if it has been done before, don't redo it, just load it
            load(data_processing_filepath);
        else % Bounds don't already exist, compute them
        %% Preprocessing (getting mean data, image bounding/ROI, filtering, fitting bounding)
            [imageData_mean,dust_filter_bounds_bottom,dust_filter_bounds_top,gate1_location_bounds,...
            gate2_location_bounds,background_totalfit,emissionlocatingdata,cutoff_height_pixels,gates,...
            delay,pixel_um_resolution,tau_fit,doublegauss_fitvariables,...
            constant_background,inclination_angle,inclination_angle_unc,nearwall_bounds,...
            near_wall_extrap,rotation_angles,total_sig_thresh,cfd_turb_prof,...
            zero_height_ref_unc,numprelim_images,imageprocesscount] = BLPreprocessing(run,...
            run_filepaths(run,:),run_start_end(run,:),run_start_end_prerun_flare(run,:),...
            runconditions_filepath,resolution_filepath,ROI(run,:),fitting_limits,synthfilepath,...
            synth_switch,synth_numimages,total_im_process,total_im_preprocess,...
            CFD_turbulent_profile_filepath,synth_real_replicate);
    
        %% Processing (primary fitting process with instantaneous velocity measurements
            [red_centroids,red_velocity,red_velocity_s,red_velocity_r,red_R2,red_g2SNR,red_signal,...
            uncertainty_wall_loc_pix,numimages] = BLImageProcessor(run,run_filepaths(run,:),run_start_end(run,:),...
            ROI(run,:),imageData_mean,dust_filter_bounds_bottom,dust_filter_bounds_top,gate1_location_bounds,...
            gate2_location_bounds,background_totalfit,emissionlocatingdata,cutoff_height_pixels,pixel_um_resolution,gates,...
            delay,fitting_limits,total_im_process,synthfilepath,synth_switch,constant_background,nearwall_bounds,...
            near_wall_extrap,rotation_angles,total_sig_thresh,cfd_turb_prof,imageprocesscount);

        %save data for next time
        save(data_processing_filepath,'red_centroids','red_velocity','red_velocity_s','red_velocity_r','red_R2','red_g2SNR',...
        'red_signal','gate1_location_bounds','gate2_location_bounds','numimages','synth_switch','inclination_angle',...
        'inclination_angle_unc','emissionlocatingdata','pixel_um_resolution','tau_fit','uncertainty_wall_loc_pix','doublegauss_fitvariables',...
        'constant_background','gate1_location_bounds', 'gate2_location_bounds','zero_height_ref_unc','background_totalfit');
    end

%% Postprocessing (filtering and calculating time-averaged velocimetry data)
[velocity_mean,velocity_mean_r,velocity_mean_s,velocity_rms,velocity_rms_r,...
    velocity_rms_s,mean_SNR,mean_R2,mean_signal,row_mm,row_mm_uncertainty,...
    sufficient_counter,filt_centroids] = BLPostprocessing(red_centroids,...
    red_velocity,red_velocity_s,red_velocity_r,red_R2,red_g2SNR,...
    red_signal,gate1_location_bounds,gate2_location_bounds,...
    numimages,synth_switch,inclination_angle(run),...
    inclination_angle_unc(run),emissionlocatingdata,pixel_um_resolution,...
    zero_height_ref_unc);

%% Plotting
BLPlotter(velocity_mean,velocity_mean_r,velocity_mean_s,velocity_rms,velocity_rms_r,...
    velocity_rms_s,mean_SNR,mean_R2,mean_signal,pixel_um_resolution,emissionlocatingdata,...
    run,numimages,tau_fit,row_mm,row_mm_uncertainty,sufficient_counter);

%% Saving
[savename] = BLSaving(velocity_mean,velocity_mean_r,...
    velocity_mean_s,velocity_rms,velocity_rms_r,...
    velocity_rms_s,mean_SNR,mean_R2,mean_signal,...
    pixel_um_resolution,emissionlocatingdata,...
    run,numimages,loc_data,row_mm,row_mm_uncertainty,...
    uncertainty_wall_loc_pix,tau_fit,...
    doublegauss_fitvariables,constant_background,...
    gate1_location_bounds,gate2_location_bounds,...
    synthfilepath,synth_switch,nondim_velo_error,...
    synth_input_tau_fit,synth_input_velocity_mean,...
    zero_height_ref_unc,sufficient_counter,...
    background_totalfit,filt_centroids);
