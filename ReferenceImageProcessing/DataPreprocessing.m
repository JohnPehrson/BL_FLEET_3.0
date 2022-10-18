clear all;close all;clc;

%% Data Preprocessing

%Clark Pehrson
%June 8, 2022

%% Load in File Path Datas
T = readtable('Data_filepaths_July22.xlsx');
runs = 1:size(T,1);
folderpaths = T{:,2};
runfolder = T{:,3};
imagefolder = T{:,4};
imagename = T{:,5};
% acefolder = T{:,8};
% acename = T{:,9};



%% Initialized
FLEET_images = zeros(length(runs),2);
x_in = zeros(2,length(runs));
y_in = zeros(2,length(runs));
g1_width = 5;
g2_width = 10;


for i = 1:length(runs)
%% Find start/end times (images) based on ACE DAQ
ace_filepath = strcat(folderpaths{i},'\',acefolder{i},'\',acename{i},'.xlsx');
DAQ_time = readmatrix(ace_filepath,'Sheet','Reduced Data','Range','A2:A1000'); %seconds
DAQ_velo = readmatrix(ace_filepath,'Sheet','Reduced Data','Range','AW2:AW1000');
    ss_times = DAQ_time(DAQ_velo>(max(DAQ_velo).*.98));
    ss_images = round(ss_times.*1000);
FLEET_images(i,:) = [min(ss_images)+3000,max(ss_images)];

%% Loop through FLEET images and find ROI data
imageData = zeros(1024,1024,2);
int_files = linspace(FLEET_images(i,1),FLEET_images(i,2),50);
for j = 1:length(int_files) %average first and last
    [imageData(:,:,j)] = Preprocessing_Imageloader(folderpaths{i},imagefolder{i},imagename{i},int_files(j));
end
imageData_avg = mean(imageData,3);

%ROI
figure;
image(imageData_avg);
colorbar;
colormap(jet(4096));
[x_loc,y_loc] = ginput(2); %top left and bottom right

x_in(:,i) = round(x_loc);
y_in(:,i) = round(y_loc);

close

%Trim ROI
ROI_imageData = imageData_avg(y_in(1,i):y_in(2,i),x_in(1,i):x_in(2,i));
ROI_imageData =  rot90(ROI_imageData,2); %rotate 180 degrees
figure;
image(ROI_imageData);
colorbar;
colormap(jet(4096));
[x_g1,y_g1] = ginput(8); %top left and bottom right
colormap(turbo(4096));
[x_g2,y_g2] = ginput(8); %top left and bottom right
close

%% Interpolate to find the general gate bounds
%g1
rows = size(ROI_imageData,1);
y_interp_1 = 0:rows;
fitobject1 = fit(y_g1,x_g1,'poly1');
x_interp_1 = fitobject1(y_interp_1);

%g2
y_interp_2 = 1:rows;
fitobject2 = fit(y_g2,x_g2,'smoothingspline');
x_interp_2 = fitobject2(y_interp_2);

x_interp_2(y_interp_2<min(y_g2)) = max(x_g2);
    %when g2 intersects g1
    gate_diff = x_interp_2-mean(x_interp_1);
    gate_cross = min(y_interp_2(gate_diff<0));
    test = y_interp_2(y_interp_2>gate_cross);

figure;
image(ROI_imageData);
colorbar;
colormap(jet(4096));
hold on;
scatter(x_g1,y_g1,'k','Linewidth',3);
hold on;
scatter(x_g2,y_g2,'r','Linewidth',3);
hold on;
plot(x_interp_1,y_interp_1,'k','Linewidth',2);
hold on;
plot(x_interp_2,y_interp_2,'r','Linewidth',2);

%% Bounds
gate1_bounds = [x_interp_1-g1_width,x_interp_1,x_interp_1+g1_width];
gate2_bounds = [x_interp_2'-g2_width,x_interp_2',x_interp_2'+g2_width];
hold on;
plot(gate1_bounds(:,1),y_interp_1,':k','Linewidth',1);
hold on;
plot(gate1_bounds(:,2),y_interp_1,':k','Linewidth',1);
hold on;
plot(gate1_bounds(:,3),y_interp_1,':k','Linewidth',1);
hold on;
plot(gate2_bounds(:,1),y_interp_2,':r','Linewidth',1);
hold on;
plot(gate2_bounds(:,2),y_interp_2,':r','Linewidth',1);
hold on;
plot(gate2_bounds(:,3),y_interp_2,':r','Linewidth',1);


end