function [gate1_location_bounds,gate2_location_bounds,time_averaged_fit,cutoff_height_pixels,...
    nearwall_bounds,background_totalfit,amplitudes,double_gauss_fitvariables,near_wall_extrap] = PrelimFitting(run,...
    imageData_mean,prerunData_mean,ypix,xpix,resolution,g1_location_col,ROI,numprelim_images,emissionlocatingdata,...
    fitting_limits,synth_switch,Delays,Gates,cfd_turb_prof,lam_run_binary,single_run,freestream_est,cross_shock_run_binary)
%this function performs preliminary fitting on the mean of the data that
%can then be applied to individual images. The primary sub-functions are to
%identify fitting bounds and observe background noise and fit it
scale = resolution/(10^6);

%% Estimate location of freestream gate 1 and 2
if cross_shock_run_binary  %downstream of the crossing shocks, can't assume freestream velocity 
    above_reflect_height = 2; %average everything above this height to estimate gate 2 location
    above_reflect_pix = emissionlocatingdata(2)-round((above_reflect_height*1000)/resolution);
    above_reflect_mean = mean(imageData_mean(1:above_reflect_pix,:),1);
    LB = [fitting_limits(1,1),0,fitting_limits(1,3),fitting_limits(1,4),g1_location_col+5,fitting_limits(1,6)];
    x0 = [fitting_limits(2,1),g1_location_col,fitting_limits(2,3),fitting_limits(2,4),g1_location_col+30,fitting_limits(2,6)];
    UB = [fitting_limits(3,1),g1_location_col+5,fitting_limits(3,3),fitting_limits(3,4),xpix,fitting_limits(3,6)];
    x = 1:xpix;
    [fitvariables] = DoubleGaussFit(x0,LB,UB,x,above_reflect_mean);
    freestream_gate_loc = [g1_location_col,fitvariables(5)]; 

else %upstream of the crossing shocks, assume roughly freestream velocity
    %find the location of the two gates in the more freestream case
    t = [ Delays(1)+Gates(1)/2 ,Delays(1)+Delays(2)+Gates(1)+Gates(2)/2];
    dt1 = (t(1))*10^(-9); %in seconds
    dx_real1 = (freestream_est.*dt1)./(scale(1));
    dt2 = (t(2)-t(1))*10^(-9); %in seconds
    dx_real2 = (freestream_est.*dt2)./(scale(1));
    dx_pix = dx_real2-dx_real1;
    freestream_gate_loc = [g1_location_col,g1_location_col+dx_pix]; 
end

%Calculating a the bounding function between the two gates as a function of
%y
if cross_shock_run_binary
    if single_run==5488 %off-axis run, very different
    cutoff_height_mm = 2; %mm, estimated height where the boundary layer becomes really substantial and emissions start to slant to zero
    else %axial runs
    cutoff_height_mm = 4; %mm, estimated height where the boundary layer becomes really substantial and emissions start to slant to zero
    end
else
    cutoff_height_mm = 8; %mm, estimated height where the boundary layer becomes really substantial and emissions start to slant to zero
end
cutoff_height_pixels = (cutoff_height_mm*1000)/resolution;
surface_height = ceil(emissionlocatingdata(2))+1; %this is the location of the wall 
cutoff_height_switchover = [round(surface_height-cutoff_height_pixels),surface_height,round(surface_height+cutoff_height_pixels)];

%Getting the bounds as a function of height using the freestream values,
%widths, and the gate locations
g1_half_width = 6;
g2_half_width = 13;
[gate1_location_bounds,gate2_location_bounds] = widthbound_fy(run,ypix,freestream_gate_loc(1),...
    freestream_gate_loc(2),cutoff_height_switchover,g1_half_width,g2_half_width,lam_run_binary);

%% Fitting Background Noise
[background_totalfit] = BackgroundFitter(run,imageData_mean,ypix,xpix,numprelim_images,synth_switch);

%plotting
close all;
linwidth = 1;
sw_dist = ((1:xpix)-freestream_gate_loc(1)).*scale*1000; %mm
vt_dist = -1*((1:ypix)-surface_height).*scale*1000; %mm
sw_dist = sw_dist';
vt_dist = vt_dist';

g1l = (gate1_location_bounds-freestream_gate_loc(1)).*scale*1000; %mm
g2l = (gate2_location_bounds-freestream_gate_loc(1)).*scale*1000; %mm

    figure(1);
    image(sw_dist,vt_dist,imageData_mean);
    set(gca,'YDir','normal');
    colorbar;
    colormap(jet(round(max(imageData_mean(:)))));
    hold on;
    plot(g1l(:,2),vt_dist,'k','Linewidth',linwidth);
    plot(g1l(:,1),vt_dist,'r','Linewidth',linwidth);
    plot(g1l(:,3),vt_dist,'r','Linewidth',linwidth);
    plot(g2l(:,1),vt_dist,'r','Linewidth',linwidth);
    plot(g2l(:,2),vt_dist,'k','Linewidth',linwidth);
    plot(g2l(:,3),vt_dist,'r','Linewidth',linwidth);
    xlabel('Streamwise Distance [mm]')
    ylabel('Vertical Distance [mm]')
    title(['Fitting Bounds for Run',num2str(single_run)])
    grid on;
    legend('Estimated Gate Location','Gate Bounds')

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

