function [noise] = NoiseCalculator(single_image_data)
    %This serves as a supplemental function to single-image fitting
    %Calculates the noise as the standard deviation of the background as
    %defined to be the lowest 50% of pixels in the image. The noise is
    %constant for the entire image.

    background_thresh = prctile(single_image_data(:),50);
    below_thresh_bin = single_image_data<background_thresh;
    background_pix = single_image_data(below_thresh_bin);
    noise = std(background_pix);

%     %% Plotting
%     figure;
%     imshow(below_thresh_bin)
%     title('Mask for Noise calculation')
%     h = gca;
%     h.Visible = 'On';
end