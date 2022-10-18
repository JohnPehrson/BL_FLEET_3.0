function [imageData_mean,dust_filter_bounds_bottom,dust_filter_bounds_top,gate1_location_bounds,...
    gate2_location_bounds,time_averaged_fit,emissionlocatingdata,cutoff_height_pixels,Gates_run,Delays_run,...
    pixel_um_resolution_run,tau_fit,doublegauss_fitvariables,...
    background_totalfit,inclination_angle,inclination_angle_unc,nearwall_bounds,...
    near_wall_extrap,rotation_angles,total_sig_thresh,cfd_turb_prof,...
    zero_height_ref_unc,numprelim_images,imageprocesscount] = BLPreprocessing(run,run_filepaths,run_timelimits,run_timelimits_flare,...
    runconditions_filepath,resolution_filepath,ROI,fitting_limits,synthfilepath,...
    synth_switch,synth_numimages,total_im_process,total_im_preprocess,CFD_turbulent_profile_filepath,...
    synth_real_replicate)
%This function will perform basic preprocessing operations related to temporal identification, 
%mean image data, image bounding/ROI, filtering, fitting bounding. Also
%calculate decay constant

%% Use data threshold
total_sig_thresh = [7500000,11000000,24000000];
total_sig_thresh_prerun = [10000000,10000000,30000000];

%% Getting Run Condition Information (resolution, gating)
load(resolution_filepath);
load(runconditions_filepath); %all time in ns

%% Image Loading and ROI if it doesn't already exist for synth and real data
if synth_switch
    numprelim_images = synth_numimages;
    imageprocesscount = round(linspace(run_timelimits(1),run_timelimits(2),numprelim_images));
    savename = ['Matfiles_preprocessing\Preprocessing_Synth_ImageData_Run',num2str(run),'_Images_',num2str(numprelim_images),'.mat'];
    if isfile(savename) %if it has been done before, don't redo it, just load it
            load(savename);
        else % Bounds don't already exist, compute them
        ROI_imagedata = []; 
        
        for imloop = imageprocesscount
            %load image
            [temp_imageData] = ImageLoader(run_filepaths(1),run_filepaths(2),synthfilepath,synth_switch,imloop);
%             %rotate 180 degrees
%             temp_ROI =  rot90(temp_ROI,2);
%             figure;
%             imagesc(temp_imageData);
            %make the temporary variable the right size for reduction
            temp_ROI_flat = temp_imageData(:);
        
            %Sort temporary variables into the reduction variables
            ROI_imagedata = [ROI_imagedata,temp_ROI_flat]; 
            
        end
        %save data for next time
        ROI_imagedata_flare = zeros(size(ROI_imagedata));
        numprelim_images_flare = total_im_preprocess;
        save(savename,'ROI_imagedata','numprelim_images_flare','ROI_imagedata_flare');
    end

else  %real data

    %% Background
        numprelim_images_flare = min([1000,;run_timelimits_flare(2)-run_timelimits_flare(1);total_im_preprocess]);
        imageprocesscount_flare = round(linspace(run_timelimits_flare(1),run_timelimits_flare(2),numprelim_images_flare));
        savename_flare = ['Matfiles_preprocessing\PreprocessingImageData_prerun_Run',num2str(run),'_Images_',num2str(numprelim_images_flare),'.mat'];
        if isfile(savename_flare) %if it has been done before, don't redo it, just load it
                load(savename_flare);
            else % Bounds don't already exist, compute them
            ROI_imagedata_flare = []; %centroid locations
            
            for imloop = imageprocesscount_flare
                %load image
                [temp_imageData] = ImageLoader(run_filepaths(1),run_filepaths(2),synthfilepath,synth_switch,imloop);
                temp_imageData = imrotate(temp_imageData,rotation_angles(run),'bilinear','crop');
%     
%                 figure;
%                 image(temp_imageData)
%                 colorbar;
%                 colormap(turbo(4096));
%                 title(['sum of intensity',num2str(round(sum(temp_imageData(:))))])

                %cut out ROI
                temp_ROI = temp_imageData(ROI(1):ROI(2),ROI(3):ROI(4));
                wall_loc_ref = zero_height_pix(1)-ROI(1);
                
%                 figure;
%                 image(temp_ROI)
%                 colorbar;
%                 colormap(turbo(4096));
                
                %make the temporary variable the right size for reduction
                temp_ROI_flat = temp_ROI(:);
            
                %Sort temporary variables into the reduction variables 
                    %if and only if the total intensity is less than some
                    %threshold value, indicating camera isn't overexposed
                    %in a lot of important places
                    if sum(temp_imageData(:))<total_sig_thresh_prerun(run)
                    ROI_imagedata_flare = [ROI_imagedata_flare,temp_ROI_flat]; 
                    end
                    numprelim_images_flare = size(ROI_imagedata_flare,2);
            end
            
            %save data for next time
            save(savename_flare,'ROI_imagedata_flare','numprelim_images_flare');
        end

    %% Real Data Information
        numprelim_images = min([1000,;run_timelimits(2)-run_timelimits(1);total_im_process]);
        image_pre_processcount = round(linspace(run_timelimits(1),run_timelimits(2),numprelim_images));
        imageprocesscount = round(linspace(run_timelimits(1),run_timelimits(2),total_im_process));
        savename = ['Matfiles_preprocessing\PreprocessingImageData_Run',num2str(run),'_Images_',num2str(total_im_process),'.mat'];
        if isfile(savename) %if it has been done before, don't redo it, just load it
                load(savename);
            else % Bounds don't already exist, compute them
            ROI_imagedata = []; %centroid locations
            ROI_imagedata_reg = [];
            
            for imloop = image_pre_processcount
                %load image
                [temp_imageData] = ImageLoader(run_filepaths(1),run_filepaths(2),synthfilepath,synth_switch,imloop);
                temp_imageData = imrotate(temp_imageData,rotation_angles(run),'bilinear','crop');
 
%                 figure;
%                 image(temp_imageData)
%                 colorbar;
%                 colormap(turbo(4096));
%                 title(['sum of intensity',num2str(round(sum(temp_imageData(:))))])
%     
                %cut out ROI
                temp_ROI = temp_imageData(ROI(1):ROI(2),ROI(3):ROI(4));
                wall_loc_ref = zero_height_pix(1)-ROI(1);
                total_ROI_int = round(sum(temp_ROI(:)));

%                 figure;
%                 image(temp_ROI)
%                 colorbar;
%                 colormap(turbo(4096));
%                 title(['sum of intensity',num2str(total_ROI_int)])

                
                %make the temporary variable the right size for reduction
                temp_ROI_flat = temp_ROI(:);
            
                %Sort temporary variables into the reduction variables 
                    %if and only if the total intensity is less than some
                    %threshold value, indicating camera isn't overexposed
                    %in a lot of important places
                    if sum(temp_ROI(:))<total_sig_thresh(run)
                    ROI_imagedata = [ROI_imagedata,temp_ROI_flat]; 
                    end
            end
            numprelim_images = size(ROI_imagedata,2);

            %save data for next time
            save(savename,'ROI_imagedata','numprelim_images');
        end
end

%% Load in the CFD turbulent boundary layer, primarily for bounding and gate-1 modifications
[cfd_turb_prof] = LoadCFD_turb_BL(CFD_turbulent_profile_filepath,pixel_um_resolution(run,1));

%% Filtering and getting the outlier bounds
totalimages = run_timelimits(1):run_timelimits(2);
[dust_filter_bounds_bottom,dust_filter_bounds_top,imageData_mean,prerunData_mean] = PixelOutlierFilter(ROI_imagedata,ROI_imagedata_flare,...
    numprelim_images,numprelim_images_flare,totalimages,0,length(ROI(1):ROI(2)),length(ROI(3):ROI(4)));

%% Getting intensity location information for time-accurate location correction
[g1_location_col,wallfit_location_col,move_flare_fit,...
    zero_height_ref_unc] = MeanLocational_Fitting(imageData_mean,length(ROI(1):ROI(2)),...
    length(ROI(3):ROI(4)),ROI,run,zero_height_pix,prerunData_mean);
emissionlocatingdata = [g1_location_col,round(wallfit_location_col),move_flare_fit];
if synth_switch
load(synth_real_replicate,'emissionlocatingdata')
end

%% Preliminaty fitting (background, local, gates) and fitting bounds
[gate1_location_bounds,gate2_location_bounds,time_averaged_fit,cutoff_height_pixels,...
    nearwall_bounds,background_totalfit,amplitudes,doublegauss_fitvariables,near_wall_extrap] = PrelimFitting(run,imageData_mean,prerunData_mean,...
    length(ROI(1):ROI(2)),length(ROI(3):ROI(4)),pixel_um_resolution(run,1),g1_location_col,ROI,...
    numprelim_images,emissionlocatingdata,fitting_limits,synth_switch,Delays,Gates,cfd_turb_prof);
background_totalfit = mean(background_totalfit(:));

%% Calculating Decay Constant
[tau_fit] = IntensityDecay(imageData_mean,time_averaged_fit,pixel_um_resolution(run,1),...
    gate1_location_bounds,gate2_location_bounds,Gates(run,:),Delays(run,:),fitting_limits,nearwall_bounds,...
    background_totalfit,amplitudes);

%% Output Variables
Gates_run = Gates(run,:);
Delays_run = Delays(run,:);
pixel_um_resolution_run = pixel_um_resolution(run,:);
end