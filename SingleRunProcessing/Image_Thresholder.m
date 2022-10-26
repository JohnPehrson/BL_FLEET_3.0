function [filtered_ROI_imagedata,thresholds] = Image_Thresholder(all_images_long,ROI,prevsrun,run)
%This function will take in a series of FLEET images and is tasked with
%generating a threshold 

sum_images_intensities = sum(all_images_long,1);
med_image_sum = median(sum_images_intensities);
lb = prctile(sum_images_intensities,2);
if prevsrun %for real data vs pre-run data, true is run data
    ub = 5*10^6+med_image_sum;
else
    ub = 2.5*10^6+med_image_sum;
    ub = min([ub,5*10^7]);
    if (run==2)||(run==3)
        lb = min(sum_images_intensities);
        ub = prctile(sum_images_intensities,8);
    end

end

%     dif_median_max = uub-med_image_sum;
%     disp(['Difference between median and 98 percentile is ',num2str(round(dif_median_max/(10^5))),'x10^5'])

thresholds = [lb,ub];
filtered_ROI_imagedata = [];
for i = 1:size(all_images_long,2)
    binary_in_bounds = (sum_images_intensities(i)>=lb)&&(sum_images_intensities(i)<=ub);
    if binary_in_bounds
        filtered_ROI_imagedata = [filtered_ROI_imagedata,all_images_long(:,i)];
    end
end

%% Stuff for plotting three 'characteristic images'
x = length(ROI(1):ROI(2));
y = length(ROI(3):ROI(4)); 

%median, at threshold, and mid-point of thresholded images
[~, min_loc] = min(abs(sum_images_intensities));
min_imdat = all_images_long(:,min_loc);
reshaped_min = reshape(min_imdat,[x,y]); 

[~, med_loc] = min(abs(sum_images_intensities-median(sum_images_intensities)));
med_imdat = all_images_long(:,med_loc);
reshaped_med = reshape(med_imdat,[x,y]); 

[~, thresh_loc] = min(abs(sum_images_intensities-ub));
thresh_imdat = all_images_long(:,thresh_loc);
reshaped_thresh = reshape(thresh_imdat,[x,y]); 

medium_filtered = mean([ub,max(sum_images_intensities)]);
[~, filt_loc] = min(abs(sum_images_intensities-medium_filtered));
filt_imdat = all_images_long(:,filt_loc);
reshaped_filt = reshape(filt_imdat,[x,y]); 


%% Visualization

figure(1);
histogram(sum_images_intensities,50);
xline(lb,'r','Linewidth',2);
xline(ub,'r','Linewidth',2);
title('Total Intensity Filtering Images')
xlabel('Total Intensity')
ylabel('Number of Images')
legend('Intensity Histogram','Limits')

figure(2);
subplot(1,4,1);
image(reshaped_min)
colorbar;
colormap(turbo(4096));
title('Min')

subplot(1,4,2);
image(reshaped_med)
colorbar;
colormap(turbo(4096));
title('Median')

subplot(1,4,3);
image(reshaped_thresh)
colorbar;
colormap(turbo(4096));
title('Max before Filtered')

subplot(1,4,4);
image(reshaped_filt)
colorbar;
colormap(turbo(4096));
title('Median of Filtered')

end