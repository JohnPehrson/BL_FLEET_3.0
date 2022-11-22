function BLPlotter(velocity_mean,velocity_mean_r,velocity_mean_s,velocity_rms,velocity_rms_r,...
    velocity_rms_s,mean_SNR,run,invar_image_compare,tau_fit,row_mm,sufficient_counter)
%This function plots the mean and root mean squared velocity using pixel
%and real scales. Optional saving of images.

%Max height
max_height = row_mm(min(sufficient_counter));

%% Setting up variables
    mean_uncertainty_combined = sqrt(velocity_mean_r.^2+velocity_mean_s.^2);
    rms_uncertainty_combined = sqrt(velocity_rms_r.^2+velocity_rms_s.^2);
    
    %% Plotting
    
    %mean and rms velocity
    figure;
    plot(velocity_mean,row_mm,'b','Linewidth',2);
    hold on;
    plot(velocity_rms,row_mm,'r','Linewidth',2);
    hold on;
    plot(velocity_mean-mean_uncertainty_combined,row_mm,'k');
    hold on;
    plot(velocity_mean+mean_uncertainty_combined,row_mm,'k');
    hold on;  
    plot(velocity_rms-rms_uncertainty_combined,row_mm,'k');
    hold on;
    plot(velocity_rms+rms_uncertainty_combined,row_mm,'k');
    hold on;
    grid on;    
    ylim([0,max_height]);
    title(['Velocimetry']);
    xlabel('Velocity [m/s]');
    ylabel('Height above the surface [mm]');
    legend('Mean Velocity','RMS Velocity');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
   
    %decay constant
    binary_rowsplot = ~isnan(velocity_mean);
    figure;
    plot(tau_fit(binary_rowsplot),row_mm(binary_rowsplot),'b','Linewidth',2);
    grid on;    
    ylim([0,max_height]);
    title('Image Decay Constant');
    xlabel('Estimated Decay Constant');
    ylabel('Height above the surface [mm]');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times

    %SNR
    figure;
    plot(mean_SNR,row_mm,'b','Linewidth',2);
    grid on;    
    ylim([0,max_height]);
    title('SNR in Gate 2');
    xlabel('SNR in Gate 2');
    ylabel('Height above the surface [mm]');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times

    %SNR vs URMS
    figure;
    scatter(mean_SNR,velocity_rms,'b','Linewidth',2);
    hold on;
    errorbar(mean_SNR,velocity_rms,rms_uncertainty_combined,'o')
    grid on;    
    title('URMS vs SNR');
    xlabel('SNR [-]');
    ylabel('URMS [m/s]');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    ylim([0,max(velocity_rms)+50])
    xlim([0,max(mean_SNR)+20])


end