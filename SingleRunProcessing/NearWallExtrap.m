function [combined_fitvars,g1_centroid_error_all,centroids_g1] = NearWallExtrap(interp_bounds,nearwall_bounds,...
    fitvariables,Delays,Gates,resolution,wall_yloc,run,snr,rows,centroid_error,cfd_turb_prof,synth_switch,near_wall_g1_scale)
%This function will identify extrapolated results in the region where
%overlapping gates is large near the surface of the test aritlce


%% Filtering the data based on the SNR to only fit good data
rows_stepper = 1:size(fitvariables,1);
snr_threshold = prctile(snr,50);            % Desired Output
binary_highsignal_not_wall = and((snr>snr_threshold),rows_stepper'>120);
binary_highsignal_not_wall(nearwall_bounds(1):nearwall_bounds(2)) = 0;

temp_fitvars = fitvariables;
y_temp = (1:rows)';
temp_fitvars = temp_fitvars(binary_highsignal_not_wall,:);
y_temp = y_temp(binary_highsignal_not_wall); %rows with good data

    %Plot amplitude and location for high snr
    high_snr_amp = temp_fitvars(:,1);
    high_snr_cent = temp_fitvars(:,2);
        all_amps = fitvariables(:,1);
        all_cents = fitvariables(:,2);
        all_rows = 1:length(all_cents);

%     figure;
%     subplot(1,2,1);
%     scatter(high_snr_amp,y_temp);
%     set(gca, 'YDir','reverse')
%     grid on;
% 
%     subplot(1,2,2);
%     scatter(high_snr_cent,y_temp);
%     set(gca, 'YDir','reverse')
%     grid on;

%% Interpolation stuff
interp_list = interp_bounds(1):interp_bounds(2);
nearwall_list = nearwall_bounds(1):nearwall_bounds(2);
above_wall_g1fit = fitvariables(interp_list,1:3);

p_loc = polyfit(y_temp,high_snr_cent,1);
p_width = polyfit(y_temp,temp_fitvars(:,3),1);
y_loc_uncorrected = polyval(p_loc,nearwall_list);
y_width = polyval(p_width,nearwall_list);

    %line estimating gate 1 above the BL
    frs_list = 1:interp_bounds(2);
    frs_wall_g1fit = polyval(p_loc,frs_list);

y_loc = y_loc_uncorrected;  %just fit as a simple line, no adjustment

%% Fitting Gate 1 Extrapolated Amplitude
    mean_amp = mean(temp_fitvars(:,1));
    all_rows = 1:max(nearwall_list);

        x0 = [mean_amp,180,30];
        LB = [100,interp_bounds(1),5];
        UB = [4000,interp_bounds(2),200];
    [fitvariables2] = SingleGaussFit(x0,LB,UB,y_temp,high_snr_amp);
    gauss_y_extrap = @(p) p(1)*exp(-(nearwall_list-p(2)).^2./(2*p(3)^2));
    gauss_y_all = @(p) p(1)*exp(-(all_rows-p(2)).^2./(2*p(3)^2));
    amp_extrap = transpose(gauss_y_extrap(fitvariables2));
    amp_all = transpose(gauss_y_all(fitvariables2));

%% Asymmetric Gauss Distribution to some pre-set amplitude
gauss_y = @(p) p(1)*exp(-(all_rows-p(2)).^2./(2*p(3)^2));
g_dist = gauss_y(fitvariables2);
cent =   round(fitvariables2(2));

%Distribution after the max
mean_max_amp = mean(amp_all(amp_all>prctile(amp_all,90))); %mean of top 90% amplitude
tail_amp = near_wall_g1_scale(1)*mean_max_amp;
diff_amp = fitvariables2(1)-tail_amp;
moved_g_dist = gauss_y([diff_amp,fitvariables2(2),fitvariables2(3).*near_wall_g1_scale(2)])+tail_amp;
g_dist_tail = g_dist;
g_dist_tail(cent:max(nearwall_list)) = moved_g_dist(cent:max(nearwall_list));

amp_extrap = g_dist_tail(nearwall_list);

        %plotting the gaussian fit
        figure;
        scatter(all_rows,all_amps,'r','Linewidth',1);
        grid on;
        hold on;
        plot(all_rows,amp_all,':c','Linewidth',2)
        plot(all_rows,g_dist_tail,'k','Linewidth',2);
        plot(nearwall_list,amp_extrap,'b','Linewidth',2)
        legend('Fit Amplitudes','Gaussian Fit','Extrapolated Amplitudes');


%% Making the output variables

    frs_wall_again = [fitvariables(frs_list,1),frs_wall_g1fit',fitvariables(frs_list,3)];
    near_wall_extrap = [amp_extrap',y_loc',y_width'];
    combined_fitvars = [frs_wall_again;near_wall_extrap];
    centroids_g1 = [combined_fitvars(:,2)];
    g1_centroid_error_all = mean(centroid_error(y_temp)).*ones(rows,1);

end