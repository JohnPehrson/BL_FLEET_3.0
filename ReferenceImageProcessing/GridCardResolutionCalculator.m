clear all;close all;clc;

%% Grid Card Resolution Calculator
%Calculates the image resolution using a mm gridcard 
%John Clark Pehrson
%March 3, 2022

%% Variables
refimage_filepath = "Data_refimages_filepaths_July22.xlsx";
T = readtable(refimage_filepath);
runs = T{:,1};
run_counter = transpose(1:length(runs));
folderpaths = T{:,2};
runfolder = T{:,3};
imagefolder = T{:,4};
imagename = T{:,5};
filepaths = T{:,6};
uniqueruns = unique(runs);
resolution_tests = strcmp(imagename,"Scale");
rotation_tests = strcmp(imagename,"Rotation");
trim_sides = [100,700,1,450;
            100,700,1,450;
            100,700,1,450];

zero_height_pix = zeros(size(resolution_tests)); %location of the surface in each run
pixel_um_resolution = zeros(length(uniqueruns),2);

%% Data wall location
mask = [300,400;300,400;300,400];
rotation_angles = [90.2,90.2,90.2];
estimate = [655;650;650];
wall_height_pix = zeros(1,length(uniqueruns));
for RUN_LOOP = 1:length(uniqueruns) % loop through different runs
        %initialize the location data for that run
            location_ims = and((runs==RUN_LOOP),rotation_tests);
            run_counter_run = run_counter(location_ims);
            wall_location_run = [];
            
        %loop through the images for each run
        for loop_uniquerun_rot = run_counter_run(1):run_counter_run(end)
        
        folderNameCat = strcat(folderpaths(loop_uniquerun_rot),"\",...
                        runfolder(loop_uniquerun_rot),"\",imagefolder(loop_uniquerun_rot)...
                        ,"\",filepaths(loop_uniquerun_rot),".tif");
        rgbImage = double(imread(folderNameCat));
        [imageData] = color_to_gray(rgbImage);
        imageData = imrotate(imageData,rotation_angles(RUN_LOOP),'bilinear','crop');
        figure;
%         subplot(1,2,1);
%         image(imageData)
%         colorbar;
%         colormap(bone(4096));
        subplot(1,2,1);
        masked_imagedata = [imageData(:,1:mask(RUN_LOOP,1)),imageData(:,mask(RUN_LOOP,2):1024)];
        image(masked_imagedata)
        colorbar;
        colormap(bone(4096));

        subplot(1,2,2);
        rowsum = sum(masked_imagedata,2);
        y = 1:length(rowsum);
        plot(rowsum,1:length(rowsum));
        set(gca, 'YDir','reverse')
        ylim([1,1024])
        hold on;

        %fitting to find the height (location of least intensity between
        %the grid card and the reflective surface)
        diff_from_est = 10;
        near_estimate_x = y((estimate(RUN_LOOP)-diff_from_est):(estimate(RUN_LOOP)+diff_from_est));
        near_estimate_y = rowsum((estimate(RUN_LOOP)-diff_from_est):(estimate(RUN_LOOP)+diff_from_est));
        [~,surf_loc] = min(near_estimate_y);
        surf_loc = estimate(RUN_LOOP)+surf_loc-1-diff_from_est;
        wall_location_run = [wall_location_run,surf_loc];
        end
    wall_height_pix(RUN_LOOP) = mean(wall_location_run);
end
close all;


%% Image resolution test
for RUN_LOOP = 1:length(uniqueruns) % loop through different runs
    location_ims = and((runs==RUN_LOOP),resolution_tests);
    run_counter_run = run_counter(location_ims);

    resolutions_run = [];
    resolutions_unc_run = [];
    for i = run_counter_run(1):run_counter_run(end)
        %% Data Loading gridcard data
        folderNameCat = strcat(folderpaths(i),"\",...
                        runfolder(i),"\",imagefolder(i)...
                        ,"\",filepaths(i),".tif");
        imageData = double(imread(folderNameCat));

        %% Trimming and averaging
        rgbImage = double(imread(folderNameCat));
        [imageData] = color_to_gray(rgbImage);
        imageData = imrotate(imageData,rotation_angles(RUN_LOOP),'bilinear','crop');
        figure;
        image(imageData)
        colorbar;
        colormap(bone(4096));

        %trim data
        imageData_trim = imageData(trim_sides(RUN_LOOP,3):trim_sides(RUN_LOOP,4),trim_sides(RUN_LOOP,1):trim_sides(RUN_LOOP,2));
        figure;
        subplot(2,2,1);
        image(imageData_trim)
        colorbar;
        colormap(bone(4096));

        %summing rows and columns
        rowsum = sum(imageData_trim,1);
        colsum = sum(imageData_trim,2);

        rowsum = rowsum./max(rowsum);
        colsum = colsum./max(colsum);

            %visualizing the row and colum intensity
            subplot(2,2,2);
            y = 1:length(colsum);
            plot(colsum,y);
            ylim([trim_sides(RUN_LOOP,3),trim_sides(RUN_LOOP,4)])
            TF_col = islocalmin(colsum,'MinProminence',0.015);
            hold on;
            plot(colsum(TF_col),y(TF_col),'r*');
            subplot(2,2,3);
            x = 1:length(rowsum);
            plot(1:length(rowsum),rowsum);
            xlim([trim_sides(RUN_LOOP,1),trim_sides(RUN_LOOP,2)])
            TF_row = islocalmin(rowsum,'MinProminence',0.015);
            hold on;
            plot(x(TF_row),rowsum(TF_row),'r*');

            gridlines_spacing_col = 1:1:length(colsum(TF_col));
            gridlines_spacing_row = 1:1:length(rowsum(TF_row));
            gridlines_col = y(TF_col);
            gridlines_row = x(TF_row);

            figure;
            scatter(gridlines_col,gridlines_spacing_col);
            hold on;
            scatter(gridlines_row,gridlines_spacing_row);
    
        %fitting to get the resolution
        p_col = polyfit(gridlines_col(2:end-1),gridlines_spacing_col(2:end-1),1);
        p_row = polyfit(gridlines_row(2:end-1),gridlines_spacing_row(2:end-1),1);
        resolutions_run = [resolutions_run,mean([p_col(1);p_row(1)])*1000]; % um/pixel

        %fitting uncertainty (resolution uncertainty)
        fitresult_col = fit(gridlines_col(2:end-1)',gridlines_spacing_col(2:end-1)','poly1');
        fitresult_row = fit(gridlines_row(2:end-1)',gridlines_spacing_row(2:end-1)','poly1');
        ci_col = confint(fitresult_col,0.95);
        ci_row = confint(fitresult_row,0.95);
        p1_uncertainty_col = ((ci_col(2,1)-ci_col(1,1))/2)*1000; % um/pixel
        p1_uncertainty_row = ((ci_row(2,1)-ci_row(1,1))/2)*1000; % um/pixel
        resolutions_unc_run = [resolutions_unc_run,mean([p1_uncertainty_col,p1_uncertainty_row])];

    end
    pixel_um_resolution(RUN_LOOP,1) = mean(resolutions_run);
    pixel_um_resolution(RUN_LOOP,2) = mean(resolutions_unc_run);

        disp(['Test ',num2str(RUN_LOOP),' has a resolution of ',num2str(pixel_um_resolution(RUN_LOOP,1)), 'um per pixel with an uncertainty of ',...
            num2str(pixel_um_resolution(RUN_LOOP,2)), 'um per pixel']);
       
        close all;
end

 %% Pointing angle downward
inclination_angle = [4.86;4.86;4.86]; %in degrees, pointing downward
inclination_angle_unc = [0.5;0.5;0.5]; %in degrees

%% Export Data
zero_height_pix = wall_height_pix;
save('C:\Users\clark\Documents\GitHub\BL_FLEET_2.0\RawDataProcessing\TestConditions/RefData.mat',...
    'pixel_um_resolution','zero_height_pix','inclination_angle','inclination_angle_unc','rotation_angles')

pixel_um_resolution(1,1)*1024/10000; %cm per length side of camera



