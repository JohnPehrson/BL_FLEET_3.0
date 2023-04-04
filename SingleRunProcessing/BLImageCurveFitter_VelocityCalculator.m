function [centroids,velocity,random_velocity_error,systematic_velocity_error,snr]= BLImageCurveFitter_VelocityCalculator(imageData_ROI,cols,...
            rows,gate1_location_bounds,gate2_location_bounds,emissionlocatingdata,fitting_limits,...
            cutoff_height_pixels,pixel_um_resolution,gates,delay,constant_background,run,cfd_turb_prof,...
            synth_switch,near_wall_g1_scale,noise)
%This function will loop through rows and apply a gaussian curve fit to it.
%The velocity is also calculated for each row. 

[centroids,snr,centroid_error,~,~,~,...
    ~,~] = Single_Image_Fitting(imageData_ROI,...
    cutoff_height_pixels,emissionlocatingdata,gate1_location_bounds,gate2_location_bounds,...
    fitting_limits,delay,gates,run,pixel_um_resolution(1),rows,cfd_turb_prof,synth_switch,...
    near_wall_g1_scale,noise);

%% Velocity Finding
velocity = zeros(rows,1);
random_velocity_error = zeros(rows,1);
systematic_velocity_error = zeros(rows,1);

for i = 1:rows
[velocity(i),random_velocity_error(i),systematic_velocity_error(i)] = VelocityFinder(centroids(i,:),centroid_error(i,:),0,pixel_um_resolution,gates,delay);
end

wall_yloc = emissionlocatingdata(2);

%velocity plot from the mean
height_plot_pix = -1.*(1:rows)+wall_yloc;
height_plot_mm = height_plot_pix.*pixel_um_resolution(1)./(10^3);

 %image viewing with bounds and the actual fit
    toppart = 1:(wall_yloc);
    figure;
    subplot(1,3,1);
    imshow(imageData_ROI(toppart,:), [], 'XData', [0 1], 'YData', [0 4]);

    subplot(1,3,2);
    plot(velocity(toppart),height_plot_mm(toppart),'k','Linewidth',2);
    ylim([0,height_plot_mm(1)]);
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlabel('Mean Velocity [m/s]');
    ylabel('Height above the surface [mm]');


    subplot(1,3,3);
    plot(snr(toppart),height_plot_mm(toppart),'k','Linewidth',2);
    ylim([0,height_plot_mm(1)]);
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlabel('Signal-to-Noise Ratio in Gate 2');
    ylabel('Height above the surface [mm]');

%% Plotting the location of centroids over the image
    figure(5);
    pix2mm = pixel_um_resolution(1)./1000;
    zeroheight = (size(imageData_ROI,1)-emissionlocatingdata(2)).*pix2mm;
    x = ((1:size(imageData_ROI,2))-gate1_location_bounds(1,2)).*pix2mm;%-3.75+137.5;
    gate_mod = -gate1_location_bounds(1,2)*pix2mm; %-3.75+137.5
    y = (size(imageData_ROI,1):-1:1).*pix2mm-zeroheight;
    image(x,y,imageData_ROI);
    axis equal;
    colormap(jet(4096));
    hold on;
    plot(centroids(:,1).*pix2mm+gate_mod,y,'--k','Linewidth',3);
    hold on;
    plot(centroids(:,2).*pix2mm+gate_mod,y,'k','Linewidth',3);
    hold on;
    plot(gate1_location_bounds(:,1).*pix2mm+gate_mod,y,'--r','Linewidth',2);
    hold on;
    plot(gate1_location_bounds(:,3).*pix2mm+gate_mod,y,'--r','Linewidth',2);
    hold on;
    plot(gate2_location_bounds(:,1).*pix2mm+gate_mod,y,'r','Linewidth',3);
    hold on;
    plot(gate2_location_bounds(:,3).*pix2mm+gate_mod,y,'r','Linewidth',3);
    axis equal;
    set(gca,'YDir','normal');
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times
    xlabel('Streamwise Distance [mm]')
    ylabel('Height above Surface [mm]')
    legend(["\mu_1","\mu_{2,m}"]);
    ylim([min(y),max(y)])
    xlim([min(x),max(x)])
    y1=ylim;        %gets current limits


%     if (run==13)
%         figure;
%         x = ((1:size(imageData_ROI,2))-gate1_location_bounds(1,2)).*pix2mm;%-3.75+137.5;
%         y = (size(imageData_ROI,1):-1:1).*pix2mm-zeroheight;
%         image(x,y,imageData_ROI);
%         axis equal;
%         colorbar;
%         colormap(bone(round(max(imageData_ROI(:)))));
%         set(gca,'YDir','normal');
%         set(gca,'FontSize', 20);
%         set(gca,'fontname','times')  % Set it to times
%         xlabel('Distance Downstream the Tripping Array [mm]')
%         ylabel('Height above the Surface [mm]')
%         ylim([min(y),max(y)])
% 
%             savename = "Scitech_Pictures/Single_Shot_Emissions.mat";
%             save(savename,'x','y','run','centroids','gate1_location_bounds','gate2_location_bounds','imageData_ROI','centroids'...
%                 ,'pix2mm','gate_mod');
% 
%     end
% 
% 
    figure(6);
    plot(snr,y','k','Linewidth',3);
    xlabel('SNR [-]')
    ylabel('Height above Surface [mm]')
    set(gca,'YDir','normal');
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times
    ylim(y1)       %sets limits to that of first subplot
    grid on;

end
