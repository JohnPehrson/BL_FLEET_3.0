function [adj,time_accurate_freestream_velo,time_av] = Syst_Decay(real,adj,synth,synth_binary,time_accurate_filepath,time_av)
%This function adjusts the mean velocity (and uncertainties) using the
%synthetic decay as a function of the decay constant.


synth_tau_corr = 50;
%% Get the velocity underprediction fit from the synthetic data (with dimensions matching the real data)
    %get synthetic data
    synth.tau_fit = synth.tau_fit+synth_tau_corr;
    fit_tau = synth.tau_fit(synth_binary);
    synth_heights = synth.velocimetry_geometricloc(synth_binary,5);
    binary_above_wall_interaction = (synth_heights>4);%&((synth_heights<11));

    %compare with the desired input velocity (from real data, but saved
    %with synthetic data)
    real_heights = real.velocimetry_geometricloc(:,5);
    real_binary = (real_heights>=min(real_heights))&(real_heights<=max(real_heights));
    real_binary(1:length(real_binary)-length(synth_heights)) = false;
    real_heights = real_heights(real_binary);
    genvelo = synth.synth_input_velocity_mean(real_binary);
    error_velo_nondim = abs(genvelo-synth.velocity_mean)./synth.velocity_mean;

    %preplotting to visualize what I should do
    figure;
    subplot(1,3,1);
    plot(genvelo,real_heights,'k','Linewidth',2);
    hold on;
    plot(synth.velocity_mean(synth_binary),real_heights,'r','Linewidth',2);
    title('Input vs. Calc Mean Velocity')
    xlabel('Mean Velocity [m/s]')
    ylabel('Height above the surface [mm]')
    legend('Input Velocity','Calculated Velocity')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    grid on;

    subplot(1,3,2);
    plot(error_velo_nondim,real_heights,'r','Linewidth',2);
    title('Error Nondimensionalized')
    xlim([0,0.1])
    xlabel('Nondimensional Error in Mean Velocity')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    grid on;

    subplot(1,3,3);
    plot(fit_tau,real_heights,'r','Linewidth',2);
    title('Decay Constant')
    xlabel('Decay Constant [ns]')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    grid on;

    %filtering for fitting
    error_velo_nondim_filt = error_velo_nondim(binary_above_wall_interaction);
    fit_tau_filt = fit_tau(binary_above_wall_interaction);

    %fitting
    nondim_velo_ftau=fit(fit_tau_filt,error_velo_nondim_filt,'poly1');
    ci = confint(nondim_velo_ftau,0.95);
    halfCI = (ci(2,:)-ci(1,:))./2;
    nondim_velo_ftau_s = @(x) sqrt((x.^2).*(halfCI(1).^2)+(halfCI(2).^2));

    %% Plotting the fit
    tau_count = linspace(min(fit_tau_filt)-10,max(fit_tau_filt)+10,200);
    figure;
    %plot the scattered points
    scatter(fit_tau_filt,error_velo_nondim_filt,'b','filled');
    hold on;
    plot(tau_count,nondim_velo_ftau(tau_count),'r','Linewidth',2);
    xlabel('Decay Constant [ns]')
    ylabel('Nondimensional Error [V_{err}/V]')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([min(fit_tau_filt)-10,max(fit_tau_filt)+10])
    title('Emission Decay Correction')

%% Using the fitting function on the real velocity
real_tau = real.tau_fit;
real_velounderprediction = nondim_velo_ftau(real_tau);
real_velounderprediction_s = nondim_velo_ftau_s(real_tau);

%% Adjust the mean velocity
adj.velocity_mean = (real_velounderprediction+1).*real.velocity_mean;
adj.velocity_mean_s = sqrt(2.*real.velocity_mean_s.^2+(real.velocity_mean.*real_velounderprediction_s).^2);

figure;
h_plot = real.velocimetry_geometricloc(:,5);
subplot(1,4,1);
plot(real.velocity_mean,h_plot);
grid on;
title('uncorrected')
subplot(1,4,2);
plot(real_tau,h_plot);
grid on;
title('tau')
subplot(1,4,3);
plot(real_velounderprediction,h_plot);
grid on;
title('correction')
subplot(1,4,4);
plot(adj.velocity_mean,h_plot);
grid on;
title('corrected')

%% Freestream mean velocity correction
load(time_accurate_filepath);
time_accurate_freestream_velo = (1+mean(real_velounderprediction)).*time_accurate_freestream_velo;

%% Adjusting time-averaged velocity
time_av_decay = time_av.decay;
% time_av_heights_bin = (time_av.height>6)&(time_av.height<12);
% time_av_decay = mean(time_av_decay(time_av_heights_bin));
time_av_uncor_mean_velo = time_av.velo;
time_av_velounderprediction = nondim_velo_ftau(time_av_decay);
time_av_cor_mean_velo = (time_av_velounderprediction+1).*time_av_uncor_mean_velo;
time_av.cor_velo = time_av_cor_mean_velo;

end