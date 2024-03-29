function [flare_data,nearwall_bounds,amplitudes,fitvariables,near_wall_extrap,...
    centroids,snr,centroid_error,y_mm,gate1_location_bounds,gate2_location_bounds,flare_height_mm] ...
    = Time_Averaged_Gate_Fitter(imageData_ROI,prerunData_mean,...
    emissionlocatingdata,cutoff_height_pixels,rows,cols,gate1_location_bounds,gate2_location_bounds,...
    fitting_limits,run,numprelim_images,synth_switch,resolution,Delays,Gates,cfd_turb_prof,sw_dist,vt_dist,...
    single_run,flare_scale,near_wall_g1_scale,create_prerun_flare_dataset,rerun_label)
%This function provides a 2d set of data that represent the first gate and
%flare data near the wall. Subtracting this 2d set of data from a
%time-resolved image should result in an image that can be fit with a
%gaussian curve near the wall easierconvertStringsToChars(rerun_label)
if synth_switch
savename = ['Preprocessing_Filestorage/NearWallFit_Synth_Run',num2str(run),'_',num2str(numprelim_images),'images.mat'];
else
savename = ['Preprocessing_Filestorage/NearWallFit_',convertStringsToChars(rerun_label),'_Run',num2str(run),'_',num2str(numprelim_images),'images.mat'];
end
    %% Check if a background fit already exists
    if isfile(savename) %if it has been done before, don't redo it, just load it
        load(savename,'flare_data','nearwall_bounds','amplitudes','fitvariables','near_wall_extrap',...
            'centroids','snr','centroid_error','y_mm','gate1_location_bounds','gate2_location_bounds','flare_height_mm');
    else % fit doesn't already exist, compute a fit here
    
    %% Prerun reference analysis to subtract time-average reflected light from the beam port
    if (~synth_switch)
    [flare_data,imageData_ROI_flaresubtracted,flare_height_mm] = Prerun_flare_processor(prerunData_mean,fitting_limits,...
        emissionlocatingdata,gate1_location_bounds,imageData_ROI,flare_scale,create_prerun_flare_dataset,...
        single_run,resolution);
    else
    flare_data = zeros(size(imageData_ROI));
    imageData_ROI_flaresubtracted = imageData_ROI;
    flare_height_pix = 0;
    flare_height_mm = 0.35;
    end


    figure(2);
    image(sw_dist,vt_dist,imageData_ROI_flaresubtracted);
    set(gca,'YDir','normal');
    colorbar;
    colormap(jet(round(max(imageData_ROI_flaresubtracted(:)))));
    xlabel('Streamwise Distance [mm]')
    ylabel('Vertical Distance [mm]')
    title(['Flare Subtraction of ',num2str(single_run)])
    grid on;

    %% Fitting the single (time-averaged) image


    %noise calculation
    [noise] = NoiseCalculator(prerunData_mean);
%     noise = max([noise,15]);

    [centroids,snr,centroid_error,fitvariables,residual,nearwall_bounds,...
        near_wall_extrap,difference_imfit] = Single_Image_Fitting(imageData_ROI_flaresubtracted,...
        cutoff_height_pixels,emissionlocatingdata,gate1_location_bounds,gate2_location_bounds,...
        fitting_limits,Delays,Gates,run,resolution,rows,cfd_turb_prof,synth_switch,near_wall_g1_scale,noise);

    %% Plotting
%     figure;
%     image(imageData_ROI_flaresubtracted)
%     colorbar;
%     colormap(jet(round(max(imageData_ROI_flaresubtracted(:)))));
%     hold on;
%     plot(gate1_location_bounds(:,1),1:rows,'r');
%     hold on
%     plot(fitvariables(:,2),1:rows,'k','Linewidth',2);
%     hold on;
%     plot(gate1_location_bounds(:,3),1:rows,'r');
%     hold on;
%     plot(gate2_location_bounds(:,1),1:rows,'r');
%     hold on
%     plot(fitvariables(:,5),1:rows,'k','Linewidth',2);
%     hold on;
%     plot(gate2_location_bounds(:,3),1:rows,'r');
%     axis equal
%     xlim([0,cols]);
%     ylim([0,400]);
%     yline(cutoff_height_pixels,'k','Linewidth',2)
%     set(gca,'FontSize', 20);
%     set(gca,'fontname','times')  % Set it to times


    %Time-averaged emissions with proper units
    x = 1:size(imageData_ROI,2);
    y = 1:size(imageData_ROI,1);
    x = x-emissionlocatingdata(1);
    y = -1.*(y-emissionlocatingdata(2));
    x_mm = resolution.*x./1000;
    y_mm = resolution.*y./1000;
    gb1 = resolution.*(gate1_location_bounds-emissionlocatingdata(1))./1000;
    gb2 = resolution.*(gate2_location_bounds-emissionlocatingdata(1))./1000;
    c1 = resolution.*(fitvariables(:,2)-emissionlocatingdata(1))./1000;
    c2 = resolution.*(fitvariables(:,5)-emissionlocatingdata(1))./1000;

    figure;
    image(x_mm,y_mm,imageData_ROI_flaresubtracted)
    colormap(jet(4096));
    hold on;
    plot(gb1(:,1),y_mm,'--r','Linewidth',2);
    plot(gb1(:,3),y_mm,'--r','Linewidth',2);
    plot(gb2(:,1),y_mm,'r','Linewidth',2);
    plot(gb2(:,3),y_mm,'r','Linewidth',2);
    plot(c1,y_mm,'--k','Linewidth',2);
    plot(c2,y_mm,'k','Linewidth',2);
    axis equal;
    set(gca, 'YDir','reverse')
    xlabel('Downstream Distance [mm]');
    ylabel('Height above the surface [mm]');
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times
    ax = gca;
    ax.YDir = 'normal';
    legend(["Fitting Bounds \mu_1","","Fitting Bounds \mu_2","","Centroids \mu_1","Centroids \mu_2"])
    xlim([min(x_mm),max(x_mm)]);
    ylim([min(y_mm),max(y_mm)]);

    figure;
    plot(snr,y_mm,'r','Linewidth',2);
    grid on;
    ylim([min(y_mm),max(y_mm)]);
    xlabel('SNR [-]');
    ylabel('Height above the surface [mm]');
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times

%Plot for visualizing near-wall fitting after extrapolation
        %     figure;
        %     image(x_mm,y_mm,difference_imfit)
        %     colorbar;
        %     colormap(jet(round(max(difference_imfit(:)))));
        %     hold on;
        %     plot(gb1(:,1),y_mm,'r','Linewidth',2);
        %     hold on
        %     plot(c1,y_mm,'k','Linewidth',2);
        %     hold on;
        %     plot(gb1(:,3),y_mm,'r','Linewidth',2);
        %     hold on;
        %     plot(gb2(:,1),y_mm,'r','Linewidth',2);
        %     hold on
        %     plot(c2,y_mm,'k','Linewidth',2);
        %     hold on;
        %     plot(gb2(:,3),y_mm,'r','Linewidth',2);
        %     axis equal;
        %     set(gca, 'YDir','reverse')
        %     xlabel('Downstream Distance [mm]');
        %     ylabel('Height above the surface [mm]');
        %     set(gca,'FontSize', 20);
        %     set(gca,'fontname','times')  % Set it to times
        %     ax = gca;
        %     ax.YDir = 'normal';
        %     legend('Emission Gate Bounds','Gate Fit')
        %     xlim([-0.5,1]);
        %     ylim([-0.25,1.25]);
        %     title('Fit overlaid on what g2 sees (extrap-subtracted)')
        
        %     figure;
        %     residual_scaled = abs(residual);
        %     image(x_mm,y_mm,residual_scaled)
        %     colorbar;
        %     colormap(jet(1000));
        %     title('Residual to the double-gaussian row-wise fit')
        %     xlabel('Downstream Distance [mm]');
        %     ylabel('Height above the surface [mm]');
        %     set(gca,'FontSize', 20);
        %     set(gca,'fontname','times')  % Set it to times
        %     ax = gca;
        %     ax.YDir = 'normal';
        %     figure;
        %     image(x_mm,y_mm,prerunData_mean)
        %     colorbar;
        %     colormap(jet(round(max(prerunData_mean(:)))));
        %     set(gca, 'YDir','reverse')
        %     xlabel('Downstream Distance [mm]');
        %     ylabel('Height above the surface [mm]');
        %     set(gca,'FontSize', 20);
        %     set(gca,'fontname','times')  % Set it to times
        %     ax = gca;
        %     ax.YDir = 'normal';
        %     title('Prerun Average')


    %% Amplitudes for the decay fitting calculation
    amplitudes = fitvariables(:,[1,4]);

  %% Saving the filter for later use if needed
      save(savename,'flare_data','nearwall_bounds','amplitudes','fitvariables','near_wall_extrap',...
          'centroids','snr','centroid_error','y_mm','gate1_location_bounds','gate2_location_bounds','flare_height_mm');
    end
end