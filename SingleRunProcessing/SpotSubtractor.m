function [filt_median,UB,LB] = SpotSubtractor(single_images_long,xlength,ylength)
%This function will identify and eliminate bright spots from instantaneous
%images. Inputs are as follows
%--median_long, an n by 1 matrix, where n is number of total pixels (x*y)
%--single_images_long, an n by m matrix, where m is number of images

%--A, represents the user tuning of level of outlier filtering, lower is
%more aggressive

%% Variable Initialization
    median_long = mean(single_images_long,2);
    total_pix = length(median_long);
    numims = size(single_images_long,2);

%% Identify outliers iteratively
    A = 2.2;               %Difference from the median the bright spots must be in terms of number of standard deviations 
    strel_size = round(2);  %pixel size to open up the filtering, accounts for the edge effects of any bright spots, must be an integer
    
    %identify rows with outliers
    kint = (single_images_long-median_long).^4;
    k_tilda = transpose((sum(kint')/(numims-1)).^0.25);
    sigint = (single_images_long-median_long).^2;
    sigma_tilda = transpose(sqrt(sum(sigint')/(numims-1)));
    binaryoutlier = k_tilda > (sigma_tilda.*A);
    indexes = find(binaryoutlier); %indexes of pixels where at least one image-pixel is an outlier
    outliercount = sum(binaryoutlier); %counts number of pixels where outliers are present
    
    while outliercount>100 %while there are still outliers, set as 100 to allow for hard-to-remove last ones

        %filter out the maximum value in the indicated rows
        for i = 1:length(indexes)
            ind = indexes(i);
            [pix_int,maxind] = max(single_images_long(ind,:));                               
            single_images_long(ind,maxind) = NaN;  %set the location of the filtered data to be NaN as a placeholder
        end

        %Recalculate the mean and re-identify rows with
        %outliers. Ignore all NaN values for now
        for i = 1:total_pix
        median_long(i) = median(single_images_long(i,~isnan(single_images_long(i,:))));          
        end
        kint = (single_images_long-median_long).^4;
        kint(isnan(kint)) = 0;  %eliminate NaN values in the summation to prevent data storage problems, but doesn't influence the summation
        k_tilda = transpose((sum(kint')/(numims-1)).^0.25);
        sigint = (single_images_long-median_long).^2;
        sigint(isnan(sigint)) = 0;  %eliminate NaN values in the summation to prevent data storage problems, but doesn't influence the summation
        sigma_tilda = transpose(sqrt(sum(sigint')/(numims-1)));
        binaryoutlier = k_tilda > (sigma_tilda.*A);
        indexes = find(binaryoutlier);
        outliercount = sum(binaryoutlier);
    end

    for i= 1:numims
        single_im_long = single_images_long(:,i);
        single_im = reshape(single_im_long,[xlength,ylength]);
        filtered_regions_im = isnan(single_im);

        %create strel object
        SE=strel('disk',strel_size);
        % apply the dilation operation.
        filtered_regions_di=imdilate(filtered_regions_im,SE);
        single_im(filtered_regions_di)=NaN;

        %put the even-more filtered single shot images back into the
        %holder-variable
        single_images_long(:,i) = single_im(:);

%         figure(1);
%         imshow(filtered_regions_di)
%         figure(2);
%         image(single_im)
%         colorbar;
%         colormap(turbo(4096));

    end

%% Replace outliers with the median
    for i = 1:total_pix
        single_images_long(i,isnan(single_images_long(i,:))) = median_long(i);
    end
    filtered_single_images_long = single_images_long;

%% Get Mean of the data, reshape it to the right size
mean_filt_image_long = mean(filtered_single_images_long,2);
filt_median = reshape(mean_filt_image_long,[xlength,ylength]);

%% Report the maximum and minimum in each pixel-location as the upper and lower bounds, respectively
    UB_long = max(filtered_single_images_long,[],2);
    UB = reshape(UB_long,[xlength,ylength]);
    LB_long = min(filtered_single_images_long,[],2);
    LB = reshape(LB_long,[xlength,ylength]);

    figure(2);
    image(filt_median);
    colorbar;
    colormap(jet(round(max(filt_median(:)))));
     
end