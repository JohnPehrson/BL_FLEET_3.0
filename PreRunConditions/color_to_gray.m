function [grayImage] = color_to_gray(rgbImage)
%This function replicates the functionality of RGB but with 12 bit images
%such as my Tif Data
 redChannel = rgbImage(:, :, 1);
              greenChannel = rgbImage(:, :, 2);
              blueChannel = rgbImage(:, :, 3);
              % Do the weighted sum.
              grayImage = .299*double(redChannel) + ...
                            .587*double(greenChannel) + ...
                            .114*double(blueChannel);
end