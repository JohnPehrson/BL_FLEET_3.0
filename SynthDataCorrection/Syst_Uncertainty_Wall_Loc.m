function [adj] = Syst_Uncertainty_Wall_Loc(real,adj,synth)
%This function propogates the uncertainty in the wall location into the
%uncertainty in the mean and RMS velocity, also in wall location

smooth_dist = 7;
h_unc = real.zero_height_ref_unc.*real.pixel_um_resolution(1).*10^(-3);
adj.velocity_mean_smooth =smooth(adj.velocity_mean,smooth_dist);
adj.velocity_mean_s_smooth =smooth(adj.velocity_mean_s,smooth_dist);
adj.velocity_rms_smooth =smooth(adj.velocity_rms,smooth_dist);
adj.velocity_rms_s_smooth =smooth(adj.velocity_rms_s,smooth_dist);
heights = real.velocimetry_geometricloc(:,5);

    for i = 2:(length(adj.velocity_mean_smooth)-1)
        height_interps = [heights(i)-h_unc,heights(i),heights(i)+h_unc];
        velos = interp1(heights,adj.velocity_mean_smooth,height_interps,'spline');
        velos_unc = interp1(heights,adj.velocity_mean_s_smooth,height_interps,'spline');
        diff_velos = [abs(velos(3)-velos(2)),abs(velos(2)-velos(1))];
        diff_unc = [abs(velos_unc(3)-velos_unc(2)),abs(velos_unc(2)-velos_unc(1))];
        add_rand_unc = max(diff_velos)+max(diff_unc);
        adj.velocity_mean_s(i) = sqrt(adj.velocity_mean_s(i).^2+add_rand_unc.^2);
    
        velos_rms = interp1(heights,adj.velocity_rms_smooth,height_interps,'spline');
        velos_rms_unc = interp1(heights,adj.velocity_rms_s_smooth,height_interps,'spline');
        diff_velos_rms = [abs(velos_rms(3)-velos_rms(2)),abs(velos_rms(2)-velos_rms(1))];
        diff_velos_rms_unc = [abs(velos_rms_unc(3)-velos_rms_unc(2)),abs(velos_rms_unc(2)-velos_rms_unc(1))];
        add_rand_unc_rms = max(diff_velos_rms)+max(diff_velos_rms_unc);
        adj.velocity_rms_s(i) = sqrt(adj.velocity_rms_s(i).^2+add_rand_unc_rms.^2);
    end


end