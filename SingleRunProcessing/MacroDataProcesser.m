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

%% Input Conditions
    %filepaths
    FLEET_folders_filepath = 'TestConditions/FLEETFilePaths.mat';          %file paths
    Run_Conditions_filepath = 'TestConditions/BLFLEETRunConditions.mat';   %stuff like gates and delays
    Resolution_filepath = 'TestConditions/RefData.mat';                    %resolution
    ACE_On_Condition_filepath = 'TestConditions/ACE_Data.mat';             %start and stops
    FLEET_ROIs_filepath = 'TestConditions/FLEET_ROIs.mat';                 %regions of interest for each run
    CFD_turbulent_profile_filepath = 'Case056NoTunnel.xlsx';

    %loading filepaths
    load(FLEET_folders_filepath);
    load(Run_Conditions_filepath);
    load(Resolution_filepath);
    load(ACE_On_Condition_filepath);
    load(FLEET_ROIs_filepath);

    %Data processing inputs
    %5479,5483,5484,5485,5486,5487,5488,5489,5490,5491,5492,5493,5494,5495
    run = [5479,5483,5484,5485,5486,5487,5488,5489,5490,5491,5492,5493,5494,5495];  %from uniqueruns, can be a vector if you want to process multiple runs in a row
    lam_runs = [5492;5493]; %these runs were closest behind the trips and should be fit with unique bounds to ensure the curve fitting algorithm can identify the gates
    cross_shock_runs = [5487;5488;5489];  %these runs were downstream the crossing shocks and thus the 'g2 bounds' should not approach the freestream
    total_im_process = NaN; %NaN if process all
    total_im_preprocess = 800; 
    lens_standoff_dist = 300; %mm
    fitting_limits =    [25,NaN,3,10,NaN,6;...
                        2000,NaN,4.5,1000,NaN,6;...
                        4096,NaN,6,4096,NaN,6]; %LB,x0,UB
    rotation_angles = rotation_angles+90; %to offset for manually rotating the images when taking scaling shots
    freestream_est_fitting = 875/(1.05); %m/s, for bounding gates
    flare_constants = [ 0,0.25,0;... 
                        0,0.15,6;...
                        0,.6,25;...
                        3,.30,0;... %5485
                        3,0.25,15;... %5486
                        4,0.15,0;...%5487
                        0,.10,5;...%5488
                        1,0.05,8;...%5489
                        3,0.15,3;...%5490
                        1,0.15,5;...%5491
                        4,0.25,0;...%5492
                        2,0.5,0;...%5493
                        0,0.05,0;...%5494
                        2,0.1,5];  %first column is moving down, second column is amplitude scale
    g1_scale_constants = [0.15,0.9;...
                        0.15,0.6;...
                        0.0,.6;...
                        .2,1;... %5485
                        0.25,.6;... %5486
                        .7,.75;...   %5487
                        0.75,.65;...   %5488
                        0.55,.65;...  %5489
                        0,.8;...   %5490
                        .5,.8;....    %5491
                        2,1;...      %5492   
                        .2,.65;...   %5493
                        0,.85;...    %5494
                        0.1,1]; %scales the amplitude of g1 near the wall. Pos for smaller G1, neg for bigger G1
    create_prerun_flare_dataset =  [true;false;false;false;false;...
                                    false;false;false;true;true;...
                                    true;true;true;true]; % create a flare-subtraction mask using the prerun data from this run?
    num_ims_runs = zeros(length(create_prerun_flare_dataset),1);

% stuff for near-wall sensitivity
    flare_amp_mod = [-0.2,-0.1,0,0.1,];
    flare_height_mod = [-1,0,1];
    flare_leftright_mod = [-4,0,4];
    g1_scale_wid = [-0.1,0,0.1];
    g1_scale_min = [-0.1,0,0.1];

    num_nw_sens = 30;
    both_const_mod = zeros(num_nw_sens,5); %[flare_constants, g1_scale_constants], list of variations
    rerun_savename = strings(num_nw_sens,1);
    for i = 1:num_nw_sens
    both_const_mod(i,:) = [randsample(flare_height_mod,1),randsample(flare_amp_mod,1),...
        randsample(flare_leftright_mod,1),randsample(g1_scale_wid,1),randsample(g1_scale_min,1)];
    rerun_savename(i) = strcat("RR",num2str(i));
    end
    rerun_switch = false;

%%General Processing Loop of images
    %specify which runs to actually use
    run_stepper = 1:length(uniqueruns);
    binary_proc_runs = false(length(uniqueruns),1);
    for i = 1:length(run)
        binary_proc_runs = or(binary_proc_runs,(run(i)==uniqueruns));
    end

for run_loop = run_stepper(binary_proc_runs)
    %% Run Details
        %getting details for that individual run
        single_run = uniqueruns(run_loop);
        run_start_end = [DAQ_start_stops(run_loop,1),DAQ_start_stops(run_loop,2)];
        run_start_end_prerun_flare = [DAQ_flare_start_stop(run_loop,1),DAQ_flare_start_stop(run_loop,2)];
        loc_data = [downstream_loc(run_loop),downstream_loc_unc(run_loop);spanwise_loc(run_loop),spanwise_loc_unc(run_loop)]; %mm, reletive the leading edge of the test article, subtract ~54 mm for this to be reletive the tripping element
        ROI = ROIs(run_loop,:);
        rot_angle = rotation_angles(run_loop);
        resolution = pixel_um_resolution(run_loop,:);
        run_inclination_angle = inclination_angle(run_loop);
        run_inclination_angle_unc = inclination_angle_unc(run_loop);
        lam_run_binary = (sum(lam_runs==uniqueruns(run_loop))==1);
        top_offset = top_offset_runs(run_loop); %for finding approx wall location
        cross_shock_run_binary = (sum(cross_shock_runs==uniqueruns(run_loop))==1);
        flare_scale = flare_constants(run_loop,:);
        near_wall_g1_scale = g1_scale_constants(run_loop,:);
        prerun_flare_binary = create_prerun_flare_dataset(run_loop);

        %total number of images to process if I want to process all of them
        if isnan(total_im_process)
            total_im_process = run_start_end(2)-run_start_end(1)-1;
        end

    %% Synthetic Data for use if the switch is flipped
        synth_switch = false;  
        synth_numimages = 500;
        
        if synth_switch
            f = fullfile('ProcessedData',['FullData_Run',num2str(run_loop),'_*']);
            fstruct = dir(f);
            synth_real_replicate = fullfile('ProcessedData',fstruct.name);
            total_im_process = synth_numimages;
    
            [synthfilepath,run_start_end_synth,nondim_velo_error,synth_input_tau_fit,...
                synth_input_velocity_mean,decay_error_eq] = Synthetic_Data_Gen(Run_Conditions_filepath,...
             Resolution_filepath,run_loop,synth_real_replicate,synth_numimages,ROI,fitting_limits);
            run_start_end = run_start_end_synth(run_loop,:);
        else
            synthfilepath = '';
            synth_numimages = total_im_process;
            nondim_velo_error = 0;
            synth_input_tau_fit = 0;
            synth_input_velocity_mean = 0;
            synth_real_replicate = 0;
            decay_error_eq = 0;
        end

    %% Image processing (if I need to do it)
        if (~synth_switch)
        data_processing_filepath = ['ProcessedData/ProcessedImageData_Run',num2str(run_loop),'Imcount_',num2str(total_im_process),'.mat'];
        else
        data_processing_filepath = ['ProcessedData/ProcessedSynthData_Run',num2str(run_loop),'Imcount_',num2str(total_im_process),'.mat'];
        end
        
%         if (isfile(data_processing_filepath)) %if it has been done before, don't redo it, just load it
%                 load(data_processing_filepath);
%                 num_ims_runs(run_loop) = numimages;
%             else % Bounds don't already exist, compute them
            %% Preprocessing (getting mean data, image bounding/ROI, filtering, fitting bounding)
                [imageData_mean,dust_filter_bounds_bottom,dust_filter_bounds_top,gate1_location_bounds,...
                gate2_location_bounds,background_totalfit,emissionlocatingdata,cutoff_height_pixels,...
                tau_fit,doublegauss_fitvariables,...
                constant_background,nearwall_bounds,...
                near_wall_extrap,mean_thresholds,prerun_thresholds,cfd_turb_prof,...
                zero_height_ref_unc,numprelim_images,imageprocesscount,flare_height_mm,y_mm] ...
                = BLPreprocessing(run_loop,...
                run_filepaths(run_loop,:),run_start_end,run_start_end_prerun_flare,...
                ROI,fitting_limits,synthfilepath,synth_switch,synth_numimages,total_im_process,total_im_preprocess,...
                CFD_turbulent_profile_filepath,synth_real_replicate,rot_angle,resolution,Gates(run_loop,:),Delays(run_loop,:),...
                lam_run_binary,single_run,top_offset,freestream_est_fitting,cross_shock_run_binary,flare_scale,near_wall_g1_scale,...
                prerun_flare_binary,...
                rerun_switch,rerun_savename,both_const_mod); %variables for rerunning the same data with different near-wall fitting parameters

        end
end
            %% Processing (primary fitting process with instantaneous velocity measurements
                [red_centroids,red_velocity,red_velocity_s,red_velocity_r,red_g2SNR,...
                uncertainty_wall_loc_pix,numimages] = BLImageProcessor(run_loop,run_filepaths(run_loop,:),...
                ROI,imageData_mean,dust_filter_bounds_bottom,dust_filter_bounds_top,gate1_location_bounds,...
                gate2_location_bounds,background_totalfit,emissionlocatingdata,cutoff_height_pixels,resolution,Gates(run_loop,:),...
                Delays(run_loop,:),fitting_limits,synthfilepath,synth_switch,constant_background,...
                rot_angle,mean_thresholds,cfd_turb_prof,imageprocesscount,near_wall_g1_scale);
    
            %save data for next time
            save(data_processing_filepath,'red_centroids','red_velocity','red_velocity_s','red_velocity_r','red_g2SNR',...
            'gate1_location_bounds','gate2_location_bounds','numimages','synth_switch','run_inclination_angle',...
            'run_inclination_angle_unc','emissionlocatingdata','pixel_um_resolution','tau_fit','uncertainty_wall_loc_pix','doublegauss_fitvariables',...
            'constant_background','gate1_location_bounds', 'gate2_location_bounds','zero_height_ref_unc','background_totalfit','y_mm','imageData_mean');
            num_ims_runs(run_loop) = numimages;


%         end

    %% Postprocessing (filtering and calculating time-averaged velocimetry data)
    [velocity_mean,velocity_mean_r,velocity_mean_s,velocity_rms,velocity_rms_r,...
        velocity_rms_s,mean_SNR,row_mm,sufficient_counter,filt_centroids,filt_velocity] ...
        = BLPostprocessing(red_centroids,...
        red_velocity,red_velocity_s,red_velocity_r,red_g2SNR,...
        gate1_location_bounds,gate2_location_bounds,...
        numimages,synth_switch,pixel_um_resolution,...
        emissionlocatingdata,y_mm,imageData_mean);
    
    %% Plotting
    BLPlotter(velocity_mean,velocity_mean_r,velocity_mean_s,velocity_rms,velocity_rms_r,...
        velocity_rms_s,mean_SNR,run_loop,numimages,tau_fit,row_mm,sufficient_counter);
    
    %% Saving
    [savename] = BLSaving(velocity_mean,velocity_mean_r,...
        velocity_mean_s,velocity_rms,velocity_rms_r,...
        velocity_rms_s,mean_SNR,...
        pixel_um_resolution,emissionlocatingdata,...
        run_loop,numimages,loc_data,row_mm,...
        uncertainty_wall_loc_pix,tau_fit,...
        doublegauss_fitvariables,constant_background,...
        gate1_location_bounds,gate2_location_bounds,...
        synthfilepath,synth_switch,nondim_velo_error,...
        synth_input_tau_fit,synth_input_velocity_mean,...
        zero_height_ref_unc,sufficient_counter,...
        background_totalfit,filt_centroids,filt_velocity,...
        decay_error_eq);
close all;

% end
