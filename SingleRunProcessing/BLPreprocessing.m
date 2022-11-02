function [imageData_mean,dust_filter_bounds_bottom,dust_filter_bounds_top,gate1_location_bounds,...
    gate2_location_bounds,time_averaged_fit,emissionlocatingdata,cutoff_height_pixels,...
    tau_fit,doublegauss_fitvariables,...
    background_totalfit,nearwall_bounds,...
    near_wall_extrap,threshs,threshs_flare,cfd_turb_prof,...
    zero_height_ref_unc,numprelim_images,imageprocesscount] = BLPreprocessing(run,run_filepaths,run_timelimits,run_timelimits_flare,...
    ROI,fitting_limits,synthfilepath,synth_switch,synth_numimages,total_im_process,total_im_preprocess,...
    CFD_turbulent_profile_filepath,synth_real_replicate,rot_angle,resolution,Gates,Delays,lam_run_binary,...
    single_run,top_offset,freestream_est,cross_shock_run_binary,flare_scale,near_wall_g1_scale,create_prerun_flare_dataset)
%This function will perform basic preprocessing operations related to temporal identification, 
%mean image data, image bounding/ROI, filtering, fitting bounding. Also
%calculate decay constant

%% Image Loading and ROI if it doesn't already exist for synth and real data
if synth_switch
    numprelim_images = synth_numimages;
    imageprocesscount = round(linspace(run_timelimits(1),run_timelimits(2),numprelim_images));
    savename = ['Preprocessing_Filestorage\Preprocessing_Synth_ImageData_Run',num2str(run),'_Images_',num2str(numprelim_images),'.mat'];
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
        savename_flare = ['Preprocessing_Filestorage\PreprocessingImageData_prerun_Run',num2str(run),'_Images_',num2str(numprelim_images_flare),'.mat'];
        if isfile(savename_flare) %if it has been done before, don't redo it, just load it
                load(savename_flare);
            else % Bounds don't already exist, compute them
            ROI_imagedata_flare = []; %centroid locations
            
            for imloop = imageprocesscount_flare
                %load image
                [temp_imageData] = ImageLoader(run_filepaths(1),run_filepaths(2),synthfilepath,synth_switch,imloop);
                temp_imageData = imrotate(temp_imageData,rot_angle,'bilinear','crop');
                
                %make the temporary variable the right size for reduction
                temp_ROI = temp_imageData(ROI(1):ROI(2),ROI(3):ROI(4));
                temp_ROI_flat = temp_ROI(:);

%                 figure(1);
%                 subplot(1,2,1);
%                 image(temp_imageData)
%                 colorbar;
%                 colormap(turbo(4096));
%                 title(['sum of intensity',num2str(round(sum(temp_imageData(:))))])
%                 subplot(1,2,2);
%                 image(temp_ROI)
%                 colorbar;
%                 colormap(turbo(4096));
%                 title(['sum of intensity',num2str(total_ROI_int)])
            
                ROI_imagedata_flare = [ROI_imagedata_flare,temp_ROI_flat]; 
            end
            [ROI_imagedata_flare,threshs_flare] = Image_Thresholder(ROI_imagedata_flare,ROI,false,run);
            numprelim_images_flare = size(ROI_imagedata_flare,2);

            %save data for next time
            save(savename_flare,'ROI_imagedata_flare','numprelim_images_flare','threshs_flare');
        end

    %% Real Data Information
        numprelim_images = min([1000,;run_timelimits(2)-run_timelimits(1);total_im_process]);
        image_pre_processcount = round(linspace(run_timelimits(1),run_timelimits(2),numprelim_images));
        imageprocesscount = round(linspace(run_timelimits(1),run_timelimits(2),total_im_process));
        savename = ['Preprocessing_Filestorage\PreprocessingImageData_Run',num2str(run),'_Images_',num2str(total_im_process),'.mat'];
        if isfile(savename) %if it has been done before, don't redo it, just load it
                load(savename);
            else % Bounds don't already exist, compute them
            ROI_imagedata = []; 
            
            for imloop = image_pre_processcount
                %load image
                [temp_imageData] = ImageLoader(run_filepaths(1),run_filepaths(2),synthfilepath,synth_switch,imloop);
                temp_imageData = imrotate(temp_imageData,rot_angle,'bilinear','crop');

                %cut out ROI
                temp_ROI = temp_imageData(ROI(1):ROI(2),ROI(3):ROI(4));
                total_ROI_int = round(sum(temp_ROI(:)));

               
%                 figure(1);
%                 subplot(1,2,1);
%                 image(temp_imageData)
%                 colorbar;
%                 colormap(turbo(4096));
%                 title(['sum of intensity',num2str(round(sum(temp_imageData(:))))])
%                 subplot(1,2,2);
%                 image(temp_ROI)
%                 colorbar;
%                 colormap(turbo(4096));
%                 title(['sum of intensity',num2str(total_ROI_int)])
                
                %make the temporary variable the right size for reduction
                temp_ROI_flat = temp_ROI(:);
                ROI_imagedata = [ROI_imagedata,temp_ROI_flat]; 

            end
            [ROI_imagedata,threshs] = Image_Thresholder(ROI_imagedata,ROI,true,run);
            numprelim_images = size(ROI_imagedata,2);

            %save data for next time
            save(savename,'ROI_imagedata','numprelim_images','threshs');
        end
end
%% Image registration with preprocessing data
xlength = length(ROI(1):ROI(2));
ylength = length(ROI(3):ROI(4));
im_reg_savename = ['Preprocessing_Filestorage\PreprocessingImageReg',num2str(run),'_RunImages_',num2str(numprelim_images),'_PreImages_',num2str(numprelim_images_flare),'.mat'];

if isfile(im_reg_savename) %if it has been done before, don't redo it, just load it
    load(im_reg_savename);
else
    [ROI_imagedata] = Preprocessing_Image_Registrar(ROI_imagedata,xlength,ylength);
    [ROI_imagedata_flare] = Preprocessing_Image_Registrar(ROI_imagedata_flare,xlength,ylength);
    save(im_reg_savename,'ROI_imagedata','ROI_imagedata_flare');
end

%% Load in the CFD turbulent boundary layer, primarily for bounding and gate-1 modifications
[cfd_turb_prof] = LoadCFD_turb_BL(CFD_turbulent_profile_filepath,resolution(1));

%% Filtering and getting the outlier bounds
savename_spots = ['Preprocessing_Filestorage\PreprocessingSpotBounds',num2str(run),'_Images_',num2str(numprelim_images_flare),'.mat'];
        if isfile(savename_spots) %if it has been done before, don't redo it, just load it
                load(savename_spots);
            else % Bounds don't already exist, compute them
                [imageData_mean,dust_filter_bounds_top,dust_filter_bounds_bottom] = SpotSubtractor(ROI_imagedata,xlength,ylength); %for run
                [prerunData_mean,~,~] = SpotSubtractor(ROI_imagedata_flare,xlength,ylength); %for prerun
                save(savename_spots,'imageData_mean','dust_filter_bounds_top','dust_filter_bounds_bottom','prerunData_mean')
        end

%% Getting intensity location information for time-accurate location correction
[g1_location_col,wallfit_location_col,move_flare_fit,zero_height_ref_unc] = ...
    MeanLocational_Fitting(imageData_mean,xlength,ylength,run,prerunData_mean,top_offset);
emissionlocatingdata = [g1_location_col,round(wallfit_location_col),move_flare_fit];
if synth_switch
load(synth_real_replicate,'emissionlocatingdata')
end

%% Preliminaty fitting (background, local, gates) and fitting bounds
[gate1_location_bounds,gate2_location_bounds,time_averaged_fit,cutoff_height_pixels,...
    nearwall_bounds,background_totalfit,amplitudes,doublegauss_fitvariables,near_wall_extrap,...
    centroids,snr,centroid_error,y_mm,gb1,gb2] = PrelimFitting(run,...
    imageData_mean,prerunData_mean,xlength,ylength,resolution(1),g1_location_col,ROI,...
    numprelim_images,emissionlocatingdata,fitting_limits,synth_switch,Delays,Gates,cfd_turb_prof,...
    lam_run_binary,single_run,freestream_est,cross_shock_run_binary,flare_scale,near_wall_g1_scale,...
    create_prerun_flare_dataset);
background_totalfit = mean(background_totalfit(:));

%% Calculating Decay Constant
[tau_fit] = IntensityDecay(imageData_mean,time_averaged_fit,resolution(1),...
    gate1_location_bounds,gate2_location_bounds,Gates,Delays,fitting_limits,nearwall_bounds,...
    background_totalfit,amplitudes);

%% Saving out time-averaged stuff
save_centroid_fit = ['ProcessedData/Time_Average_Fit_Run',num2str(run),'_',num2str(numprelim_images),'images.mat'];
  save(save_centroid_fit,'run','centroids','snr','centroid_error','numprelim_images','y_mm','gb1','gb2','tau_fit','zero_height_ref_unc');

end