clear all;close all;clc;

%% Grid Card Resolution Calculator
%Calculates the image resolution using a mm gridcard 
%John Clark Pehrson
%March 3, 2022

%% Variables
refimage_filepath = "RefImages_Filepaths.xlsx";
T = readtable(refimage_filepath);
runs = T{:,1};
run_counter = transpose(1:length(runs));
folderpaths = T{:,2};
runfolder = T{:,3};
imagefolder = T{:,4};
imagetype = T{:,5};
imagename = T{:,6};
imagestarts = T{:,7};
imageends = T{:,8};
uniqueruns = unique(runs);
resolution_tests = strcmp(imagetype,"Scale");
rotation_tests = strcmp(imagetype,"Rotation");

%Per-run information that isn't easily standardized because the camera was
%moving between runs
scale_ROI = [200,525,500,1024;
            100,500,500,1024;
            100,525,425,960;
            100,525,400,1024;
            100,500,300,1024;
            150,480,1,525;
            100,480,200,800;
            150,475,1,575;
            150,450,125,700;
            100,420,175,875;
            75,460,260,1024;
            150,470,280,1024;
            150,460,60,625;
            125,440,375,950];

rot_ROI = [700,900,400,1024;
           600,800,400,1024;
           700,850,400,1024;
           700,850,400,1024;
           650,850,300,1024;
           650,800,150,700;
           600,800,100,900;
           650,800,100,700;
           600,800,1,700;
           600,750,1,800;
           650,800,1,1024;
           650,800,1,1024;
           600,800,1,450;
           600,800,350,900];

%Initializing Variables for the runs
wall_height_pix = zeros(1,length(uniqueruns));
image_rotation_degrees = zeros(1,length(uniqueruns));
wall_surfs = zeros(length(uniqueruns),2);
pixel_um_resolution = zeros(length(uniqueruns),2);

%% Data wall location and rotation
for RUN_LOOP = 1:length(uniqueruns) % loop through different runs
        %initialize the location data for that run
            location_ims = and((runs==uniqueruns(RUN_LOOP)),rotation_tests);
            image_numbers_run = imagestarts(location_ims):imageends(location_ims);
            
        %variables to track the rotation and wall location for duplicate
        %ref images per tunnel run    
        run_rotations = [];
        run_wall_locations = [];
        run_wall_surf = [];

        %loop through the images for each run
        for loop_uniquerun_rot = image_numbers_run
        
        folderNameCat = strcat(folderpaths(location_ims),"\",...
                        runfolder(location_ims),"\",imagefolder(location_ims)...
                        ,"\",imagename(location_ims),num2str(loop_uniquerun_rot),"\",...
                        imagename(location_ims),num2str(loop_uniquerun_rot),".tif");

        %Calculations
            %load in and trim ref image
            imageData = double(imread(folderNameCat));
            masked_imagedata = imageData(rot_ROI(RUN_LOOP,1):rot_ROI(RUN_LOOP,2),rot_ROI(RUN_LOOP,3):rot_ROI(RUN_LOOP,4));
            %row-averaging to visualize the surface
            rowsum = sum(masked_imagedata,2);
            y_plot = 1:length(rowsum);
            %near-wall gradient of intensity
            [Gmag,Gdir] = imgradient(masked_imagedata);
            %turning the image into a binary map
            P_thresh = prctile(Gmag(:),96.5);
            BW = imbinarize(Gmag,P_thresh);
            [y,x] = find(BW);  %// get x y coordinates of all curve points
            %find the lower half of the binary points, which should represent
            %the gradient of the light at the surface, which is the 'wall'
                y_thresh = prctile(y,50);
                lower = y>y_thresh;
                xl = x(lower);
                yl = y(lower);
                xu = x(~lower);
                yu = y(~lower);
                BW(yu,xu) = 0;
                [y,x] = find(BW);  %// get x y coordinates of what I expect to be the near-wall gradient
                ymedian = median(y);
                    y_bound = 10;
                yrange = [ymedian-y_bound,ymedian+y_bound];
                passed_points = (y>yrange(1))&(y<yrange(2));
                xp = x(passed_points);
                yp = y(passed_points);
                xf = x(~passed_points);
                yf = y(~passed_points);
                BW(yf,xf) = 0;

            wallsurffit = polyfit( xp, yp, 1 );  %// fit 3rd deg poly
            xx = 1:size(BW,2);
            wallfit_yy = polyval( wallsurffit, xx );
            %moving the gradient down 1 pixel towards where the wall actually
            %is (as the gradient is more of the taper-off than the 'wall')
            ROI_mod_wallfit_yy = wallfit_yy+1;
            %modying to sit in the entire image from the ROI
            xx_full = 1:1024;
            wallsurffit_full = [wallsurffit(1),wallsurffit(2)+1+rot_ROI(RUN_LOOP,1)-1];
            yy_full =polyval(wallsurffit_full,xx_full);

%         %Plotting
%         figure(1);
%         subplot(1,3,1);
%         image(imageData)
%         colorbar;
%         colormap(bone(max(imageData(:))));
%         hold on;
%         plot(xx_full, yy_full, '.-r', 'LineWidth', 2 );
%         hold off;
% 
%         subplot(1,3,2);
%         image(masked_imagedata)
%         colorbar;
%         colormap(bone(max(imageData(:))));
%         title(strcat(runfolder(location_ims)," image number ",num2str(loop_uniquerun_rot)))
%         hold on;
%         plot(xx, ROI_mod_wallfit_yy, '.-r', 'LineWidth', 2 );
%         hold off;
% 
%         subplot(1,3,3);
%         plot(rowsum,y_plot);
%         set(gca, 'YDir','reverse')
%         ylim([1,max(y_plot)])
% 
%         %Image Gradient
%         figure(2);
%         subplot(1,2,1);
%         image(Gmag)
%         colorbar;
%         hold on;
%         plot(xx, wallfit_yy, '.-r', 'LineWidth', 2 );
%         hold off;
%    
%         subplot(1,2,2);
%         imshow(BW, 'border', 'tight' );
%         hold all
%         plot(xx, wallfit_yy, '.-', 'LineWidth', 2 );
%         h = gca;
%         h.Visible = 'On';

        run_rotations = [run_rotations,atand(wallsurffit(1))];
        run_wall_locations= [run_wall_locations,mean(yy_full)];
        run_wall_surf = [run_wall_surf;wallsurffit_full];
        end
    image_rotation_degrees(RUN_LOOP) = mean(run_rotations);
    wall_height_pix(RUN_LOOP) = mean(run_wall_locations);
    wall_surfs(RUN_LOOP,:) = mean(run_wall_surf,1);
end
close all;


%% Image resolution test
for RUN_LOOP = 1:length(uniqueruns) % loop through different runs
    location_ims = and((runs==uniqueruns(RUN_LOOP)),resolution_tests);
    image_numbers_run = imagestarts(location_ims):imageends(location_ims);

    resolutions_run = [];
    resolutions_unc_run = [];

    for loop_uniquerun_rot = image_numbers_run
        %% Data Loading gridcard data
        folderNameCat = strcat(folderpaths(location_ims),"\",...
                        runfolder(location_ims),"\",imagefolder(location_ims)...
                        ,"\",imagename(location_ims),num2str(loop_uniquerun_rot),"\",...
                        imagename(location_ims),num2str(loop_uniquerun_rot),".tif");

        %% Trimming and averaging
        imageData = double(imread(folderNameCat));
        imageData = imrotate(imageData,image_rotation_degrees(RUN_LOOP),'bilinear','crop');

        %trim data
        imageData_trim = imageData(scale_ROI(RUN_LOOP,1):scale_ROI(RUN_LOOP,2),scale_ROI(RUN_LOOP,3):scale_ROI(RUN_LOOP,4));

        %summing rows and columns
        rowsum = sum(imageData_trim,1);
        colsum = sum(imageData_trim,2);

        rowsum = rowsum./max(rowsum);
        colsum = colsum./max(colsum);

        TF_row = islocalmin(rowsum,'MinProminence',0.015);
        TF_col = islocalmin(colsum,'MinProminence',0.015);
        gridlines_spacing_col = 1:1:length(colsum(TF_col));
        gridlines_spacing_row = 1:1:length(rowsum(TF_row));
            x = 1:length(rowsum);
            y = 1:length(colsum);
        gridlines_col = y(TF_col);
        gridlines_row = x(TF_row);
    
%         %Plotting
%             figure(1);
%             image(imageData)
%             colorbar;
%             colormap(bone(round(max(imageData(:)))));
%     
%             figure(2);
%             subplot(2,2,1);
%             hold off;
%             image(imageData_trim)
%             colormap(bone(round(max(imageData(:)))));
% 
%             %visualizing the row and colum intensity
%             subplot(2,2,2);
%             hold off;
%             plot(colsum,y);
%             ylim([1,length(colsum)])
%             hold on;
%             plot(colsum(TF_col),y(TF_col),'r*');
%             set(gca, 'YDir','reverse')
%             title('Column Sum')
% 
%             subplot(2,2,3);
%             hold off;
%             plot(1:length(rowsum),rowsum);
%             xlim([1,length(rowsum)])
%             hold on;
%             plot(x(TF_row),rowsum(TF_row),'r*');
%             title('Row Sum')
% 
%             subplot(2,2,4);
%             hold off;
%             scatter(gridlines_col,gridlines_spacing_col);
%             hold on;
%             scatter(gridlines_row,gridlines_spacing_row);
%             legend('Column Sum','Row Sum')
%     
        %fitting to get the resolution
        p_col = polyfit(gridlines_col(2:end-1),gridlines_spacing_col(2:end-1),1);
        p_row = polyfit(gridlines_row(2:end-1),gridlines_spacing_row(2:end-1),1);
        if length(gridlines_col)>length(gridlines_row)
            resolutions_run = [resolutions_run,p_col(1)*1000]; % um/pixel
                    
            %fitting uncertainty (resolution uncertainty)
            fitresult_col = fit(gridlines_col(2:end-1)',gridlines_spacing_col(2:end-1)','poly1');
            ci_col = confint(fitresult_col,0.95);
            p1_uncertainty_col = ((ci_col(2,1)-ci_col(1,1))/2)*1000; % um/pixel
            resolutions_unc_run = [resolutions_unc_run,p1_uncertainty_col,p1_uncertainty_row];

        else
            resolutions_run = [resolutions_run,p_row(1)*1000]; % um/pixel

            %fitting uncertainty (resolution uncertainty)
            fitresult_row = fit(gridlines_row(2:end-1)',gridlines_spacing_row(2:end-1)','poly1');
            ci_row = confint(fitresult_row,0.95);
            p1_uncertainty_row = ((ci_row(2,1)-ci_row(1,1))/2)*1000; % um/pixel
            resolutions_unc_run = [resolutions_unc_run,p1_uncertainty_row];

        end

    end
    pixel_um_resolution(RUN_LOOP,1) = mean(resolutions_run);

    if std(resolutions_run)>1
        disp(['At least one run doesnt match up well in ',num2str(uniqueruns(RUN_LOOP))]);
    end

    pixel_um_resolution(RUN_LOOP,2) = mean(resolutions_unc_run);

        disp(['Run ',num2str(uniqueruns(RUN_LOOP)),' has a resolution of ',num2str(pixel_um_resolution(RUN_LOOP,1)), 'um per pixel with an uncertainty of ',...
            num2str(pixel_um_resolution(RUN_LOOP,2)), 'um per pixel']);

        close all;
end

 %% Pointing angle downward
inclination_angle = 2.*ones(1,length(uniqueruns)); %in degrees, pointing downward
inclination_angle_unc = 0.1.*ones(1,length(uniqueruns)); %in degrees

%% Export Data
zero_height_pix = wall_height_pix;
rotation_angles = image_rotation_degrees;
save('C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/RefData.mat',...
    'uniqueruns','pixel_um_resolution','zero_height_pix','inclination_angle','inclination_angle_unc',...
    'rotation_angles','wall_surfs')

pixel_um_resolution(1,1)*1024/10000; %cm per length side of camera



