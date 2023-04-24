clear all;close all;clc;

%% FLEET ROI Picker
    %This script loops through each set of FLEET images, pulls an average,
    %and then lets me manually set in points in the rough area I want to
    %curve fit

 %% Get info about directories to know where to load in data
   currentdir  = pwd;
 idcs   = strfind(currentdir,'\');
 rootdir = currentdir(1:idcs(end)-1);
savefilepath = fullfile(rootdir,"SingleRunProcessing","TestConditions");

%% Input variables from other preprocessing scripts

    FLEET_folders_filepath = fullfile(savefilepath,"FLEETFilePaths.mat");
    Run_Conditions_filepath = fullfile(savefilepath,"BLFLEETRunConditions.mat");    %stuff like gates and delays
    ACE_On_Condition_filepath = fullfile(savefilepath,"ACE_Data.mat");
    Resolution_filepath = fullfile(savefilepath,"RefData.mat");

    %loading filepaths
    load(FLEET_folders_filepath);
    load(Run_Conditions_filepath);
    load(ACE_On_Condition_filepath);
    load(Resolution_filepath);

    %total number of images to load in from each run
    tot_run_im = 100;
    num_runs = size(run_filepaths,1);
    pixs_side = 1024;

    %initialization
    single_shot_images = zeros(pixs_side^2,tot_run_im);
    averaged_images = zeros(pixs_side,pixs_side,num_runs);
    ROIs = zeros(num_runs,4);
    est_wall_loc = zeros(num_runs,1);

%% Loop through each image and load, average
all_fleet_averages = 'FLEET_averages_100.mat';
if isfile(all_fleet_averages) %if it has been done before, don't redo it, just load it
    load(all_fleet_averages);
else
    
    for i = 1:num_runs
        folderName = run_filepaths(i,1);
        imageName = run_filepaths(i,2);
        imageNumbers = round(linspace(DAQ_start_stops(i,1),DAQ_start_stops(i,2),tot_run_im));
        rot_angle = rotation_angles(i)+90;
    
        for j = 1:tot_run_im
            image_filepath = strcat(folderName,'\',imageName,num2str(imageNumbers(j),'%06.f'),".tif");
            single_image_data_untrim = double(imread(image_filepath));
            single_shot_images(:,j) = single_image_data_untrim(:);
    
        end
            long_average = median(single_shot_images,2);
            reshaped_average = reshape(long_average,[pixs_side,pixs_side]); 
            reshaped_average_rot = imrotate(reshaped_average,rot_angle,'bilinear','crop');
            averaged_images(:,:,i) = reshaped_average_rot;
    end

save(all_fleet_averages,'averaged_images');    
end

%% Loop through averaged FLEET images, manually prescribe the location of emissions
FLEET_gate_locations = 'FLEET_gate_locations.mat';
g1_ims = 4;
g2_ims = 8;
x_g1_s = zeros(num_runs,g1_ims);
y_g1_s = zeros(num_runs,g1_ims);
x_g2_s = zeros(num_runs,g2_ims);
y_g2_s = zeros(num_runs,g2_ims);

if isfile(FLEET_gate_locations) %if it has been done before, don't redo it, just load it
    load(FLEET_gate_locations);
else
    for i = 1:num_runs
        FLEET_average_image = averaged_images(:,:,i);
        
        %picking points
        figure(1);
        image(FLEET_average_image)
        colormap(jet(round(max(FLEET_average_image(:)))));
        disp('Pick 4 points on the first gate')
        [x_g1,y_g1] = ginput(g1_ims);
        disp('Pick 8 points on the second gate')
        [x_g2,y_g2] = ginput(g2_ims);
    
        %vis points
        figure(1);
        image(FLEET_average_image)
        hold on;
        plot(x_g1,y_g1,'r','Linewidth',2);
        plot(x_g2,y_g2,'r','Linewidth',2)

        %save points 
        x_g1_s(i,:) = x_g1;
        y_g1_s(i,:) = y_g1;
        x_g2_s(i,:) = x_g2;
        y_g2_s(i,:) = y_g2;

    end
    
    save(FLEET_gate_locations,'x_g1_s','y_g1_s','x_g2_s','y_g2_s')
end

%Picking a ROI based on the manually picked gate locations
top_offset_runs = zeros(1,num_runs);
for i = 1:num_runs
    x_g1 = x_g1_s(i,:);
    y_g1 = y_g1_s(i,:);
    x_g2 = x_g2_s(i,:);
    y_g2 = y_g2_s(i,:);

    bottom_offset = 25;
    front_offset = 50;
    back_offset = 50;
    if (i~=1)&&(i<10)
        top_offset = 300;
    elseif i==1
        top_offset = 500;
    else
        top_offset = 350;
    end
    top_offset_runs(i) = top_offset;

    %lower bound lower than the lowest point
    minheight = max([y_g1,y_g2]);
    lower = minheight+bottom_offset;
    upper = minheight-top_offset;
    left = mean(x_g1)-front_offset;
    right = max([x_g1,x_g2])+back_offset;

    ROIs(i,:) = round([upper,lower,left,right]); %y then x
end

%visualizing everything put together
for i = 1:num_runs 
    FLEET_average_image = averaged_images(:,:,i);
    ROI = ROIs(i,:);
    FLEET_average_trimmed = FLEET_average_image(ROI(1):ROI(2),ROI(3):ROI(4));

    figure(1);
    image(FLEET_average_trimmed)
    colormap(jet(round(max(FLEET_average_trimmed(:)))));
 
end

%% Save out final ROIs
ROI_savefile = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/FLEET_ROIs.mat';
save(ROI_savefile,'ROIs','top_offset_runs');


