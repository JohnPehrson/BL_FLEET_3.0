function [combined_fitvars,g1_centroid_error_all,centroids_g1] = NearWallExtrap(interp_bounds,nearwall_bounds,...
    fitvariables,Delays,Gates,resolution,wall_yloc,run,snr,rows,centroid_error,cfd_turb_prof,synth_switch)
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
y_temp = y_temp(binary_highsignal_not_wall);

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


%% Correct y_loc using the velocity profile and known stuff about the flow

v_inf = 868;
time_int = (Delays(1)+Gates(1)/2).*(10^-9);
max_pix_disp = v_inf*time_int/(resolution*10^-6);

if (run==1)  %laminar, estimate with piecewise linear
    scale_to_wall = (nearwall_list-wall_yloc)./min((nearwall_list-wall_yloc));
    scale_to_wall(scale_to_wall<0) = 0;
    scale_to_wall = (scale_to_wall-1).*-1;
    y_loc = y_loc_uncorrected-scale_to_wall*max_pix_disp;
else %better approximated with a turbulent boundary layer
    nearwall_list_adjusted = -1.*(nearwall_list-wall_yloc);
    nearwall_list_adjusted(nearwall_list_adjusted<0) = 0;

    scale_to_wall = interp1(cfd_turb_prof.y_pix,cfd_turb_prof.v_nondim,nearwall_list_adjusted,'spline');
    scale_to_wall = (scale_to_wall-1).*-1;
    if (run==2)
        scale_to_wall = scale_to_wall.*.1;
    elseif (run==3)
        scale_to_wall = scale_to_wall.*0.4;
    end
    y_loc = y_loc_uncorrected-scale_to_wall*max_pix_disp;

end


% figure;
% plot(cfd_turb_prof.v_nondim,cfd_turb_prof.y_pix)

    %% gaussian amplitude fitting         
    amp_g1 = above_wall_g1fit(:,1);
    amp_g1_sorted = sort(amp_g1,'descend');
    if run==3
        amp_scale = 1;
        mean_amp = amp_g1_sorted(round(length(amp_g1_sorted)/amp_scale));
        mean_amp = mean_amp.*0.9;
    else
        amp_scale = 1.2;
        mean_amp = amp_g1_sorted(round(length(amp_g1_sorted)/amp_scale));
    end

    if synth_switch
        mean_amp = mean_amp.*0.55;
    end


%                 x0 = [2500,180,30];
%                 LB = [100,interp_bounds(1),5];
%                 UB = [4000,interp_bounds(2),200];
%             [fitvariables] = SingleGaussFit(x0,LB,UB,interp_list',amp_g1);
%             gauss_y = @(p) p(1)*exp(-(nearwall_list-p(2)).^2./(2*p(3)^2));
%             amp_extrap = gauss_y(fitvariables);

    frs_wall_again = [fitvariables(frs_list,1),frs_wall_g1fit',fitvariables(frs_list,3)];
    near_wall_extrap = [mean_amp.*ones(size(y_loc')),y_loc',y_width'];
    combined_fitvars = [frs_wall_again;near_wall_extrap];
    centroids_g1 = [combined_fitvars(:,2)];

%% Extrapolate the uncertainty in the linear fit of the centroid for the first gate
g1_centroid_error_all = mean(centroid_error(y_temp)).*ones(rows,1);


end