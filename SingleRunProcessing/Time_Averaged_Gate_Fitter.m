function [flare_data,nearwall_bounds,amplitudes,fitvariables,near_wall_extrap] = Time_Averaged_Gate_Fitter(imageData_ROI,prerunData_mean,...
    emissionlocatingdata,cutoff_height_pixels,rows,cols,gate1_location_bounds,gate2_location_bounds,...
    fitting_limits,run,numprelim_images,synth_switch,resolution,background_totalfit,Delays,Gates,cfd_turb_prof)
%This function provides a 2d set of data that represent the first gate and
%flare data near the wall. Subtracting this 2d set of data from a
%time-resolved image should result in an image that can be fit with a
%gaussian curve near the wall easier
if synth_switch
savename = ['Matfiles_preprocessing/NearWallFit_Synth_Run',num2str(run),'_',num2str(numprelim_images),'images.mat'];
else
savename = ['Matfiles_preprocessing/NearWallFit_Run',num2str(run),'_',num2str(numprelim_images),'images.mat'];
end
    %% Check if a background fit already exists
    if isfile(savename) %if it has been done before, don't redo it, just load it
        load(savename,'flare_data','nearwall_bounds','amplitudes','fitvariables','near_wall_extrap');
    else % fit doesn't already exist, compute a fit here
    
    %% Prerun reference analysis to subtract time-average reflected light from the beam port
    if (~synth_switch)
    [flare_data,imageData_ROI_flaresubtracted] = Prerun_flare_processor(prerunData_mean,fitting_limits,emissionlocatingdata,gate1_location_bounds,run,imageData_ROI);
    else
    flare_data = zeros(size(imageData_ROI));
    imageData_ROI_flaresubtracted = imageData_ROI;
    end

    %% Fitting the single (time-averaged) image
    [centroids,snr,centroid_error,R2,fitvariables,residual,nearwall_bounds,...
        near_wall_extrap] = Single_Image_Fitting(imageData_ROI_flaresubtracted,...
        cutoff_height_pixels,emissionlocatingdata,gate1_location_bounds,gate2_location_bounds,...
        fitting_limits,Delays(run,:),Gates(run,:),run,resolution,rows,cfd_turb_prof,synth_switch);

    %% Plotting

titles_runs = ["Laminar","Pizza Box","Synthetic Roughness"];

figure; 
image(imageData_ROI_flaresubtracted);
colorbar;
colormap(bone(4096));
set(gca, 'YDir','reverse')
hold on;
plot(centroids(:,1),1:rows,':r','Linewidth',3);
plot(centroids(:,2),1:rows,'r','Linewidth',3);
set(gca,'FontSize', 20);
set(gca,'fontname','times')  % Set it to times
legend('Initial Emissions Centroids','Displaced Emissions Centroids');

    figure;
    image(imageData_ROI_flaresubtracted)
    colorbar;
    colormap(jet(round(max(imageData_ROI_flaresubtracted(:)))));
    hold on;
    plot(gate1_location_bounds(:,1),1:rows,'r');
    hold on
    plot(fitvariables(:,2),1:rows,'b');
    hold on;
    plot(gate1_location_bounds(:,3),1:rows,'r');
    hold on;
    plot(gate2_location_bounds(:,1),1:rows,'r');
    hold on
    plot(fitvariables(:,5),1:rows,'b');
    hold on;
    plot(gate2_location_bounds(:,3),1:rows,'r');
    axis equal
    xlim([0,cols]);
    ylim([0,400]);
    title(titles_runs(run))

    figure;
    residual_scaled = (residual./max(abs(residual(:)))+1).*124;
    image(residual_scaled)
    colorbar;
    colormap(jet(256));
    title('Residual to the double-gaussian row-wise fit')

    %Time-averaged emissions with proper units
    x = 1:size(imageData_ROI,2);
    y = 1:size(imageData_ROI,1);
    x = x-emissionlocatingdata(1);
    y = -1.*(y-emissionlocatingdata(2));
    x_mm = resolution.*x./1000;
    y_mm = resolution.*y./1000;

    figure;
    image(x_mm,y_mm,imageData_ROI)
    colorbar;
    colormap(turbo(1500));
    axis equal;
    set(gca, 'YDir','reverse')
        if synth_switch
        title(['Synthetic Data']);
        elseif run ==1
        title(['Laminar']);
        elseif run ==2
        title(['Pizza Box']);
        else %run = 3
        title(['Synthetic Roughness']);
        end
    xlabel('Downstream Distance [mm]');
    ylabel('Height above the surface [mm]');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    ax = gca;
    ax.YDir = 'normal';

%% Amplitudes for the decay fitting calculation
amplitudes = fitvariables(:,[1,4]);

  %% Saving the filter for later use if needed
      save(savename,'flare_data','nearwall_bounds','amplitudes','fitvariables','near_wall_extrap');
    end
end