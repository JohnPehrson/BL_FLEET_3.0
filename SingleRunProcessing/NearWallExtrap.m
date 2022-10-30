function [combined_fitvars,g1_centroid_error_all,centroids_g1] = NearWallExtrap(interp_bounds,nearwall_bounds,...
    fitvariables,Delays,Gates,resolution,wall_yloc,run,snr,rows,centroid_error,cfd_turb_prof,synth_switch,near_wall_g1_scale)
%This function will identify extrapolated results in the region where
%overlapping gates is large near the surface of the test aritlce


%% Filtering the data based on the SNR to only fit good data
snr_sort = sort(snr(interp_bounds(1):interp_bounds(2),1),'descend');                                  % Sort Descending for the SNR of gate 1
snr_threshold = min(snr_sort(1:ceil(length(snr_sort)*0.75)));             % Desired Output
binary_highsignal_not_wall = (snr>snr_threshold);
binary_highsignal_not_wall(nearwall_bounds(1):nearwall_bounds(2)) = 0;

temp_fitvars = fitvariables;
y_temp = (1:rows)';
temp_fitvars = temp_fitvars(binary_highsignal_not_wall,:);
y_temp = y_temp(binary_highsignal_not_wall); %rows with good data

%% Interpolation stuff
interp_list = interp_bounds(1):interp_bounds(2);
nearwall_list = nearwall_bounds(1):nearwall_bounds(2);
above_wall_g1fit = fitvariables(interp_list,1:3);

p_loc = polyfit(y_temp,temp_fitvars(:,2),1);
p_width = polyfit(y_temp,temp_fitvars(:,3),1);
y_loc_uncorrected = polyval(p_loc,nearwall_list);
y_width = polyval(p_width,nearwall_list);

    %line estimating gate 1 above the BL
    frs_list = 1:interp_bounds(2);
    frs_wall_g1fit = polyval(p_loc,frs_list);

y_loc = y_loc_uncorrected;  %just fit as a simple line, no adjustment

%% Fitting Gate 1 Extrapolated Amplitude
    mean_amp = mean(temp_fitvars(:,1));
    
        x0 = [mean_amp,180,30];
        LB = [100,interp_bounds(1),5];
        UB = [4000,interp_bounds(2),200];
        amps = temp_fitvars(:,1);
    [fitvariables2] = SingleGaussFit(x0,LB,UB,y_temp,amps);
    gauss_y = @(p) p(1)*exp(-(nearwall_list-p(2)).^2./(2*p(3)^2));
    amp_extrap = transpose(gauss_y(fitvariables2));
    
    %Scaling by a manual factor
    dist_from_wall = transpose(abs(wall_yloc-nearwall_list));
    linear_scaling = dist_from_wall/max(dist_from_wall);
    max_lin = max(linear_scaling);
    exp_scaling = abs(linear_scaling).^2;
    mult_fact = (1-near_wall_g1_scale)+exp_scaling*near_wall_g1_scale;

    frs_wall_again = [fitvariables(frs_list,1),frs_wall_g1fit',fitvariables(frs_list,3)];
    near_wall_extrap = [amp_extrap.*mult_fact,y_loc',y_width'];
    combined_fitvars = [frs_wall_again;near_wall_extrap];
    centroids_g1 = [combined_fitvars(:,2)];


%% Extrapolate the uncertainty in the linear fit of the centroid for the first gate
g1_centroid_error_all = mean(centroid_error(y_temp)).*ones(rows,1);


end