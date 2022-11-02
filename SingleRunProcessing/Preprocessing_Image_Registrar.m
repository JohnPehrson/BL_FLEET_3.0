function [reg_time_ac_images] = Preprocessing_Image_Registrar(time_ac_images,xlength,ylength)
%This function will register the images reletive their mean to try and
%provide a more precise time-accurate mean intensity

%% Variables
num_ims = size(time_ac_images,2);
reg_time_ac_images = zeros(size(time_ac_images));
move_im = zeros(num_ims,2);
%% Get the unfiltered, unregistered image intensity mean
mean_time_ac_images = median(time_ac_images,2);
mean_time_s = reshape(mean_time_ac_images,[xlength,ylength]);

% figure;
% image(mean_time_s);
% colorbar;
% colormap(jet(round(max(mean_time_s(:)))));
% title('Unregistered Mean')
% 

point_loc = [60,250];

%% Loop through all images, register them
num_ims = size(time_ac_images,2);
for i= 1:num_ims
    single_im = time_ac_images(:,i);
    single_im_s = reshape(single_im,[xlength,ylength]);
    [single_im_s_reg,move_im(i,:)] = ImageRegistration(single_im_s,mean_time_s,point_loc);
    reg_time_ac_images(:,i) = single_im_s_reg(:);
end

%% Get the mean again, compare against the original (is this worth doing at all?)
% reg_mean_time_ac_images = median(reg_time_ac_images,2);
% reg_mean_time_s = reshape(reg_mean_time_ac_images,[xlength,ylength]);
% 
% figure;
% image(reg_mean_time_s);
% colorbar;
% colormap(jet(round(max(reg_mean_time_s(:)))));
% title('Registered Mean')
% 
% figure;
% imshowpair(mean_time_s,reg_mean_time_s,"diff")
% colorbar;
% colormap(jet);
end