clear all;close all;clc;

%% Video Maker (Main)
% This script turns a series of processed images into a video file

run = 2;
folderpath = "C:\Users\clark\Documents\GitHub\BL_FLEET_2.0\RawDataProcessing\Fit_images\R";
real_filename = "Tot_im";
filt_filename = "filtcents_ims";
bounds = [100,400];
num_im_override = 100;

%% Get the folder with the most images
[real_filepath] = ImageFilepath(folderpath,run,real_filename);

%% Load in the images
[image_data_all,centroids_all,numims] = LoadImages(real_filepath);

if ~isnan(num_im_override) %if override
numims = num_im_override;
end


%% Load in the filtered centroids
[filt_centroids] = Load_filtered_centroids(folderpath,run,filt_filename);

%% Video Creater
y = 1:size(centroids_all,1);
filt_binary = ~isnan(filt_centroids);

vid_filepath = strcat("Vid_nofit_Run",num2str(run),"_Ims",num2str(numims),".avi");
v = VideoWriter(vid_filepath,'Uncompressed AVI');
v.FrameRate = 20;
open(v)
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
        y_c2_bin = (y_c2>bounds(1))&(y_c2<bounds(2));
        y_c2 = y_c2(y_c2_bin)-bounds(1)+1;
        c2 = c2(y_c2_bin);

    image(image_matrix(bounds(1):bounds(2),:));
    colormap(bone(2048));
    colorbar;
    axis equal;
    set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
    
%     hold on;
%     plot(c1(bounds(1):bounds(2)),y_c1,'r','Linewidth',2);
%     scatter(c2,y_c2,'*r','Linewidth',2);
%     hold off;

    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v)