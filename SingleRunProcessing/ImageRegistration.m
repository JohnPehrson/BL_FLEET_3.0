function [imageData_ROI_translated,move_ROI_ud] = ImageRegistration(imageData_ROI,imageData_mean,emissionlocatingdata)

%This function will find the row associated with the wall location and then
%move the image to match the wall location with some time-consistent wall
%location

%% Image Registration
[optimizer, metric] = imregconfig('monomodal');
imageData_ROI_translated = imregister(imageData_ROI,imageData_mean,'translation',optimizer,metric);

%% Tracking the translations
tform = imregtform(imageData_ROI,imageData_mean,'translation',optimizer,metric);

x = emissionlocatingdata(1);
y = emissionlocatingdata(2);
[x_trans,y_trans] = transformPointsForward(tform,x,y);
x_rel = x_trans-x;
y_rel = y_trans-y;

%% Move image ROI
%imageData_ROI_translated
move_ROI_ud = [x_rel,y_rel];

end