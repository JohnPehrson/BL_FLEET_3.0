function [gate1_location_bounds,gate2_location_bounds,time_averaged_fit,cutoff_height_pixels,...
    nearwall_bounds,background_totalfit,amplitudes,double_gauss_fitvariables,near_wall_extrap] = PrelimFitting(run,...
    imageData_mean,prerunData_mean,ypix,xpix,resolution,g1_location_col,ROI,numprelim_images,emissionlocatingdata,...
    fitting_limits,synth_switch,Delays,Gates,cfd_turb_prof)
%this function performs preliminary fitting on the mean of the data that
%can then be applied to individual images. The primary sub-functions are to
%identify fitting bounds and observe background noise and fit it

%% Finding the location of gate 1
rowaveraged_data = mean(imageData_mean);
x = 1:xpix;
x0 = [1000,g1_location_col,4];
LB = [100,25,1];
UB = [4096,xpix-25,8];
[gate1_fitvariables] = SingleGaussFit(x0,LB,UB,x,rowaveraged_data);
freestream_gate_loc = [gate1_fitvariables(2)];

%% Calculating gaussian fitting bounds
%find the location of the two gates in the more freestream case
TF = islocalmax(imageData_mean,2,'MinProminence',25);
numlocalmax = sum(TF);
x0 = [ypix/2,gate1_fitvariables(2),2,ypix/5,gate1_fitvariables(2)+25,4];
LB = [10,25,1,10,gate1_fitvariables(2),1];
UB = [ypix,xpix-25,4,ypix,xpix-25,8];
[gate2_fitvariables] = DoubleGaussFit(x0,LB,UB,x,numlocalmax);
freestream_gate_loc = [freestream_gate_loc,gate2_fitvariables(5)+2]; %the plus 2 is to very approximately account for decay

%Calculating a the bounding function between the two gates as a function of
%y
cutoff_height_mm = [3,3.5,4]; %mm
cutoff_height_mm = cutoff_height_mm(run); %mm, estimated height where the boundary layer becomes really substantial and emissions start to slant to zero
cutoff_height_pixels = (cutoff_height_mm*1000)/resolution;
surface_height = ceil(emissionlocatingdata(2))+1; %this is the location of the wall 
cutoff_height_switchover = [round(surface_height-cutoff_height_pixels),surface_height,round(surface_height+cutoff_height_pixels)];

%Getting the bounds as a function of height using the freestream values,
%widths, and the gate locations
g1_half_width = 5;
g2_half_width = 9;
[gate1_location_bounds,gate2_location_bounds] = widthbound_fy(run,ypix,freestream_gate_loc(1),freestream_gate_loc(2),cutoff_height_switchover,g1_half_width,g2_half_width);

%% Fitting Background Noise
[background_totalfit] = BackgroundFitter(run,imageData_mean,ypix,xpix,numprelim_images,synth_switch);

%% Fitting the flare and providing a subtraction data set for near wall velocimetry
imageData_mean_nobackground = imageData_mean-background_totalfit;
[flare_g1_fit,nearwall_bounds,amplitudes,double_gauss_fitvariables,near_wall_extrap] = Time_Averaged_Gate_Fitter(imageData_mean_nobackground,...
prerunData_mean,emissionlocatingdata,cutoff_height_pixels,...
ypix,xpix,gate1_location_bounds,gate2_location_bounds,fitting_limits,run,numprelim_images,synth_switch,...
resolution,background_totalfit,Delays,Gates,cfd_turb_prof);
time_averaged_fit = background_totalfit+flare_g1_fit;

    %plotting
    figure;
    image(imageData_mean-time_averaged_fit)
    colorbar;
    colormap(bone(4096));
    hold on;
    plot(gate1_location_bounds(:,1),1:ypix,'r');
    hold on
    plot(gate1_location_bounds(:,2),1:ypix,'b');
    hold on;
    plot(gate1_location_bounds(:,3),1:ypix,'r');
    hold on;
    plot(gate2_location_bounds(:,1),1:ypix,'r');
    hold on
    plot(gate2_location_bounds(:,2),1:ypix,'b');
    hold on;
    plot(gate2_location_bounds(:,3),1:ypix,'r');

%     %% Making a single shot vs time-averaged figure for journal paper 
%     %load single shot
%     single_shot_filepath = "C:\Users\clark\Documents\Grad_Research\Data\Feb_FLEET\Test1_Feb_C001H001S0001\Test1_Feb_C001H001S0001019853.tif";
%     singleshot_image = double(imread(single_shot_filepath));
%     singleshot_ROI = singleshot_image(ROI(1):ROI(2),ROI(3):ROI(4));
%     singleshot_ROI =  rot90(singleshot_ROI,2);
%     
%     %trim data to only be above the surface
%     imageData_mean_subtracted = imageData_mean-time_averaged_fit;
%     singleshot_subtracted = singleshot_ROI-time_averaged_fit;
%     imageData_mean_abovesurface = imageData_mean;
%     singleshot_ROI_abovesurface = singleshot_ROI;
%     imageData_mean_subtracted_abovesurface = imageData_mean_subtracted;
%     singleshot_subtracted = singleshot_subtracted;
%     
%     pix2mm = resolution/1000;
%     zeroheight = (size(singleshot_ROI_abovesurface,1)-(emissionlocatingdata(2))+45).*pix2mm;
%     x = (1:size(singleshot_ROI_abovesurface,2)).*pix2mm+135;
%     y = (size(singleshot_ROI_abovesurface,1):-1:1).*pix2mm-zeroheight;
%     
%     figure;
%     subplot(1,2,1);
%     image(x,y,singleshot_subtracted)
%     set(gca,'YDir','normal');
%     colormap(hot(4096));
%     set(gca,'FontSize', 20);
%     set(gca,'fontname','times')  % Set it to times
% 
%     subplot(1,2,2);
%     image(x,y,imageData_mean_subtracted_abovesurface)
%     set(gca,'YDir','normal');
%     colormap(hot(4096));
%     set(gca,'FontSize', 20);
%     set(gca,'fontname','times')  % Set it to times
%     xlabel('Distance Downstream the Tripping Array [mm]')
% 
% 
%  
%     figure;
%     image(x,y,imageData_mean_abovesurface)
%     set(gca,'YDir','normal');
%     colormap(hot(4096));
%     colorbar;
%     set(gca,'FontSize', 20);
%     set(gca,'fontname','times')  % Set it to times
%     xlabel('Distance Downstream the Tripping Array [mm]')
%     ylabel('Height above the Surface [mm]')

% 
% %Time averaged fitting to look for image artefacts
% rows = 90:160;
% sumrows_mean = mean(imageData_mean(rows,:));
% figure;
% plot(1:length(sumrows_mean),sumrows_mean);
% 
end

