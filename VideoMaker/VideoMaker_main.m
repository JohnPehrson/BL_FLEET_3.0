clear all;close all;clc;

%% Video Maker (Main)
% This script turns a series of processed images into a video file

run = 8;
folderpath = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\Fit_images\R";
centroid_folderpath =  "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\Filtered_ProcessedData";
wall_loc_folderpath = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\ProcessedData";
real_filename = "Tot_im";
centroid_filename = strcat("FilteredImageData_Run",num2str(run));
wall_loc_data_filename = "FullData_Run";
res_filepath = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions\RefData.mat";
load(res_filepath);

num_im_override = 500;

%% Get the folder with the most images
[real_filepath] = ImageFilepath(folderpath,run,real_filename);

%% Get the wall locations from the full, filtered dataset
[y_mm,rows_plot] = Load_FullDataset(wall_loc_folderpath,run,wall_loc_data_filename);


%% Load in the images
[image_data_all,centroids_all,numims] = LoadImages(real_filepath);

if ~isnan(num_im_override) %if override
numims = num_im_override;
end

%% Load in the filtered centroids
[filt_centroids] = Load_filtered_centroids(centroid_folderpath,run,centroid_filename);
bounds = [1,size(filt_centroids,1)];

    %extrapolate all pixel wall locations
    all_pix = 1:size(filt_centroids,1);
    p = polyfit(rows_plot,y_mm,1);
    y_mm_all = polyval(p,all_pix);

    x_all = 1:size(image_data_all,2);
    x_mm_all = (pixel_um_resolution(run,1)./1000).*x_all;

%% Video Creater
y = 1:size(centroids_all,1);
filt_binary = ~isnan(filt_centroids);

vid_filepath = strcat("Vid_nofit_Run",num2str(run),"_Ims",num2str(numims),".avi");
v = VideoWriter(vid_filepath,'Uncompressed AVI');
v.FrameRate = 20;
% open(v)
fullfig;
for i = 1:numims

    image_matrix = image_data_all(:,:,i);
    centroids = centroids_all(:,:,i);
        c1 = centroids(:,1);
        c2 = centroids(:,2);
        filt_centroids = filt_binary(:,(2*i-1):2*i);
        c2_bin = filt_centroids(:,2);
        c2 = c2(c2_bin);
        y_c2 = y(c2_bin);
        y_c1 = y(bounds(1):bounds(2))-bounds(1)+1;
        y_c2_bin = (y_c2>=bounds(1))&(y_c2<=bounds(2));
        y_c2 = y_c2(y_c2_bin)-bounds(1)+1;
        c2 = c2(y_c2_bin);
        y_mm_all_c2 = y_mm_all(c2_bin);

    image(x_mm_all,y_mm_all,image_matrix);
    colormap(bone(4096));
    colorbar;
    axis equal;
    xlabel('Pixels');
    ylabel('Height above the surface [mm]');
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times
    ax = gca;
    ax.YDir = 'normal';
    xlim([min(x_mm_all),max(x_mm_all)]);
    ylim([min(y_mm),max(y_mm)]);
    
    hold on;
    c1_mm = (pixel_um_resolution(run,1)./1000).*c1;
    c2_mm = (pixel_um_resolution(run,1)./1000).*c2;
    plot(c1_mm(bounds(1):bounds(2)),y_mm_all,'r','Linewidth',2);
    scatter(c2_mm,y_mm_all_c2','*r','Linewidth',2);
    hold off;

    pause(0.5)
%     frame = getframe(gcf);
%     writeVideo(v,frame);
end

% close(v)