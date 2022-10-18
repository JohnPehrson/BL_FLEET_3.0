close all;clear all;clc;

%% Synthetic data correction program
%John Clark Pehrson
%March 18, 2022

%This script compares the real and synthetic replication of perpendicular
%FLEET testing and then applies some corrections to the real data

%% Initializing Variables
run_titles = ["Laminar";"Pizza Box";"Synthetic Roughness"];
lens_standoff_dist = 363; %mm
center_mag_height = 8; %mm
resolution_filepath = 'C:\Users\clark\Documents\GitHub\BL_FLEET_2.0\RawDataProcessing\TestConditions\RefData.mat';%resolution

run = [1,2,3];
run_i = 2;
%  for i = run
    folderpath = 'C:\Users\clark\Documents\GitHub\BL_FLEET_2.0\RawDataProcessing\ProcessedData';
    real_filename = strcat("ProcessedData_Run",num2str(run(run_i)),"_Images");
    real_stringlength = strlength(real_filename);
    synth_filename = strcat("ProcessedData_Synth_Run",num2str(run(run_i)),"_Images");
    synth_stringlength = strlength(synth_filename);
    time_accurate_filepath = "TimeAccurateFLEET\timeaccuratedata.mat";
    
    %laminar data case
    lam_real_filename = strcat("ProcessedData_Run",num2str(1),"_Images");
    lam_real_stringlength = strlength(real_filename);
    lam_synth_filename = strcat("ProcessedData_Synth_Run",num2str(1),"_Images");
    lam_synth_stringlength = strlength(synth_filename);

    %time-averaged stuff
    preproc_filepath = "C:\Users\clark\Documents\GitHub\BL_FLEET_2.0\RawDataProcessing\Matfiles_preprocessing";
    centroids_filepath = strcat("NearWallFit_Run",num2str(run(run_i)),"_");
    runconditions_filepath = 'C:\Users\clark\Documents\GitHub\BL_FLEET_2.0\RawDataProcessing\TestConditions\BLFLEETRunConditions.mat';    %stuff like gates and delays
    
    %% Data Loading
    [real,synth,maxrealimages] = RealSynthDataLoader(folderpath,...
        real_filename,real_stringlength,synth_filename,synth_stringlength);
    
    %% Time-averaged Data Loading
    [time_av] = TimeAv_Centroid_Loading(preproc_filepath,centroids_filepath,...
        real.pixel_um_resolution,runconditions_filepath,run(run_i),real.emissionlocatingdata);
    
    %% Making a new data set represnting data after systematic adjustment
    adj = real;
    realheights_binary = real.velocimetry_geometricloc(:,5)>=0;
    realheights = real.velocimetry_geometricloc(realheights_binary,5);
    synthheights_binary = (synth.velocimetry_geometricloc(:,5)>=0);
    synthheights = synth.velocimetry_geometricloc(synthheights_binary,5);
    
    %% Mean Velocity Systematic Adjustment
    [adj,time_accurate_freestream_velo,time_av] = Syst_Decay(real,adj,synth,synthheights_binary,time_accurate_filepath,time_av);

    %% RMS Velocity Systematic Adjustment
    [adj] = Syst_Precision(real,synth,adj,synthheights_binary,realheights_binary,...
        folderpath,lam_real_filename,lam_real_stringlength,lam_synth_filename,...
        lam_synth_stringlength,run(run_i));
%  end
close all;

%% Propogating wall location uncertainty into mean and RMS velocity (widening uncertainty bounds near the wall)
[adj] = Syst_Uncertainty_Wall_Loc(real,adj,synth);

%% Correcting for the systematic error in mean/rms velocity due to inclination (magnification of the line unequally)
%find the difference in distance based on inclination
L_hvec = adj.velocimetry_geometricloc(:,5)-center_mag_height;
load(resolution_filepath);
angle_beam = 3; %likely beam angle, degrees
d_vec = L_hvec.*sind(inclination_angle(run_i)+angle_beam);
Mo_Md = 1-d_vec./lens_standoff_dist;
adj.velocity_mean = adj.velocity_mean.*Mo_Md;
adj.velocity_rms = adj.velocity_rms.*Mo_Md;


L_hvec_av = time_av.height-center_mag_height;
d_vec_av = L_hvec_av.*sind(inclination_angle(run_i)+angle_beam);
Mo_Md = 1-d_vec_av./lens_standoff_dist;
time_av.cor_velo = time_av.cor_velo.*Mo_Md;

%% Raw data plotting vs synthetic data
    real.mean_uncertainty_combined = sqrt(real.velocity_mean_r.^2+real.velocity_mean_s.^2);
    real.rms_uncertainty_combined = sqrt(real.velocity_rms_r.^2+real.velocity_rms_s.^2);
    synth.mean_uncertainty_combined = sqrt(synth.velocity_mean_r.^2+synth.velocity_mean_s.^2);
    synth.rms_uncertainty_combined = sqrt(synth.velocity_rms_r.^2+synth.velocity_rms_s.^2);
    max_h = max(realheights);

    %mean and rms velocity
    figure;
    plot(real.velocity_mean(realheights_binary),realheights,'b','Linewidth',2);
    hold on;
    plot(real.velocity_rms(realheights_binary),realheights,'r','Linewidth',2);
    hold on;
    plot(time_av.velo,time_av.height,'k','Linewidth',2);
    hold on;
    plot(real.velocity_mean(realheights_binary)-real.mean_uncertainty_combined(realheights_binary),realheights,'k');
    hold on;
    plot(real.velocity_mean(realheights_binary)+real.mean_uncertainty_combined(realheights_binary),realheights,'k');
    hold on;  
    plot(real.velocity_rms(realheights_binary)-real.rms_uncertainty_combined(realheights_binary),realheights,'k');
    hold on;
    plot(real.velocity_rms(realheights_binary)+real.rms_uncertainty_combined(realheights_binary),realheights,'k');
    hold on;
    grid on;    
    ylim([0,max_h]);
    title(['Velocimetry before Synthetic Correction']);
    xlabel('Velocity [m/s]');
    ylabel('Height above the surface [mm]');
    legend('Mean Velocity','RMS Velocity','Time-Averaged Mean Velo');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times

    hold on;
    plot(synth.velocity_mean(synthheights_binary),synthheights,'g','Linewidth',2); 
    hold on;
    plot(synth.velocity_rms(synthheights_binary),synthheights,'m','Linewidth',2);
    hold on;
    plot(synth.velocity_mean(synthheights_binary)-synth.mean_uncertainty_combined(synthheights_binary),synthheights,'k');
    hold on;
    plot(synth.velocity_mean(synthheights_binary)+synth.mean_uncertainty_combined(synthheights_binary),synthheights,'k');
    hold on;  
    plot(synth.velocity_rms(synthheights_binary)-synth.rms_uncertainty_combined(synthheights_binary),synthheights,'k');
    hold on;
    plot(synth.velocity_rms(synthheights_binary)+synth.rms_uncertainty_combined(synthheights_binary),synthheights,'k');
    hold on;
    grid on;    
    ylim([0,max_h]);
    title(['Velocimetry before Synthetic Correction']);
    xlabel('Velocity [m/s]');
    ylabel('Height above the surface [mm]');
    legend('Mean Velocity','RMS Velocity');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times

    %mean SNR and R2 
    figure;
    plot(real.mean_SNR(realheights_binary),realheights,'b','Linewidth',2);
    hold on;
    plot(real.mean_R2(realheights_binary),realheights,'r','Linewidth',2);
    grid on;    
    ylim([0,max_h]);
    title('Mean SNR and R2');
    xlabel('[-]');
    ylabel('Height above the surface [mm]');
    legend('Mean SNR','Mean R2');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times

    hold on;
    plot(synth.mean_SNR(synthheights_binary),synthheights,'g','Linewidth',2);
    hold on;
    plot(synth.mean_R2(synthheights_binary),synthheights,'m','Linewidth',2);
    grid on;    
    ylim([0,max_h]);
    title('Mean SNR and R2');
    xlabel('[-]');
    ylabel('Height above the surface [mm]');
    legend('Mean SNR','Mean R2');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times

    %decay constant
    figure;
    plot(real.tau_fit(realheights_binary),realheights,'b','Linewidth',2);
    grid on;    
    ylim([0,max_h]);
    title('Image Decay Constant');
    xlabel('Estimated Decay Constant');
    ylabel('Height above the surface [mm]');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times

    hold on;
    plot(synth.tau_fit(synthheights_binary),synthheights,'g','Linewidth',2);
    grid on;    
    ylim([0,max_h]);
    title('Image Decay Constant');
    xlabel('Estimated Decay Constant');
    ylabel('Height above the surface [mm]');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times

%% Plotting adjusted data (after synthetic data corrections)
adj_binary = adj.velocimetry_geometricloc(:,5)>=0;
adjheights = adj.velocimetry_geometricloc(:,5);
adj.mean_uncertainty_combined = sqrt(adj.velocity_mean_r.^2+adj.velocity_mean_s.^2);
adj.rms_uncertainty_combined = sqrt(adj.velocity_rms_r.^2+adj.velocity_rms_s.^2);

%mean correction
figure;
plot(real.velocity_mean(realheights_binary),realheights,'r','Linewidth',2.5);
hold on;
plot(real.velocity_mean(realheights_binary)-real.mean_uncertainty_combined(realheights_binary),realheights,':r','Linewidth',1.5);
hold on; 
plot(real.velocity_mean(realheights_binary)+real.mean_uncertainty_combined(realheights_binary),realheights,':r','Linewidth',1.5);
hold on;
plot(adj.velocity_mean_smooth,adjheights,'k','Linewidth',2.5);
hold on;
plot(adj.velocity_mean_smooth-adj.mean_uncertainty_combined,adjheights,':k','Linewidth',1.5);
hold on; 
plot(adj.velocity_mean_smooth+adj.mean_uncertainty_combined,adjheights,':k','Linewidth',1.5);
hold on;  
grid on;
xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
ylabel('Height above the surface [mm]');
legend('Uncorrected','95% CI','95% CI','Decay-corrected','95% CI','95% CI');
set(gca,'FontSize', 20);
set(gca,'fontname','times')  % Set it to times
xlim([0,1000]);
title('Emission Decay Correction')

%% Side-by-side Mean and RMS
figure;
plot(adj.velocity_mean_smooth,adjheights,'k','Linewidth',2.5);
hold on;
plot(adj.velocity_mean_smooth-adj.mean_uncertainty_combined,adjheights,':k','Linewidth',1.5);
hold on; 
%uncorrected
plot(real.velocity_mean(realheights_binary),realheights,'r','Linewidth',1.5);
hold on;
%timeav
plot(time_av.cor_velo,time_av.height,'b','Linewidth',1.5);
hold on;
%unc. cont.
plot(adj.velocity_mean_smooth+adj.mean_uncertainty_combined,adjheights,':k','Linewidth',1.5);
hold on;  
grid on;
xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
ylabel('Height above the surface [mm]');
legend('Corrected','95% CI','Uncorrected','Time-average Mean Velocity');
set(gca,'FontSize', 20);
set(gca,'fontname','times')  % Set it to times
xlim([0,1000]);
ylim([0,max(adjheights)])
title(run_titles(run_i))

figure;
frs_binary = adjheights>10;
frstrm = mean(adj.velocity_mean_smooth(frs_binary));
adj.velocity_rms_perc = adj.velocity_rms_smooth./frstrm;
adj.rms_uncertainty_combined_perc = adj.rms_uncertainty_combined./frstrm;

plot(adj.velocity_rms_perc,adjheights,'k','Linewidth',2.5);
hold on;
plot(adj.velocity_rms_perc-adj.rms_uncertainty_combined_perc,adjheights,':k','Linewidth',1.5);
hold on;
plot(real.velocity_rms(realheights_binary)./frstrm,realheights,'r','Linewidth',1.5);
hold on;
plot(real.velocity_rms(realheights_binary)./frstrm+real.rms_uncertainty_combined(realheights_binary)./frstrm,realheights,':r','Linewidth',1.5);
hold on;
plot(real.velocity_rms(realheights_binary)./frstrm-real.rms_uncertainty_combined(realheights_binary)./frstrm,realheights,':r','Linewidth',1.5);
hold on;
plot(adj.velocity_rms_perc+adj.rms_uncertainty_combined_perc,adjheights,':k','Linewidth',1.5);
grid on;    
xlabel('$V_{RMS} / \bar{V}_\infty$ [m/s]','Interpreter','Latex')
ylabel('Height above the surface [mm]');
legend('Corrected','95% CI','Uncorrected','95% CI');
set(gca,'FontSize', 20);
set(gca,'fontname','times')  % Set it to times
title(run_titles(run_i))

figure;
adj.mean_SNR_smooth = smooth(adj.mean_SNR,5);
plot(adj.mean_SNR_smooth,adjheights,'k','Linewidth',2.5);
xlabel('SNR');
ylabel('Height above the surface [mm]');
set(gca,'FontSize', 20);
set(gca,'fontname','times')  % Set it to times

%% Plotting adjusted data and comparing with CFD
%load in cfd
cfdfilepath = "C:\Users\clark\Documents\GitHub\BL_FLEET_2.0\SyntheticDataCorrection\Case063-5-1.txt";
cfd = readtable(cfdfilepath);
cfd_y = cfd{:,3}.*1000;
cfd_velo = cfd{:,4};

    %mean velocity with CFD (RANS)
    figure;
    plot(adj.velocity_mean_smooth,adjheights,'k','Linewidth',2.5);
    hold on;
    plot(time_av.cor_velo,time_av.height,'g','Linewidth',2);
    hold on;
    plot(adj.velocity_mean_smooth-adj.mean_uncertainty_combined,adjheights,':k','Linewidth',1.5);
    hold on;
    plot(cfd_velo,cfd_y,'r','Linewidth',2);
    hold on;
    plot(adj.velocity_mean_smooth+adj.mean_uncertainty_combined,adjheights,':k','Linewidth',1.5);
    grid on;    
    ylim([0,10]);
    xlim([0,1000]);
    xlabel('Velocity [m/s]');
    ylabel('Height above the surface [mm]');
    legend('Time-accurate Mean Velocity','Time-average Mean Velocity',...
        '95% Confidence Interval','RANS CFD Mean Velocity');
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times

    %% Plotting comparison with RMS
    %load in cfd
    cfdfilepath = "C:\Users\clark\Documents\GitHub\BL_FLEET\SyntheticDataCorrection\Case063-5-1.txt";
    cfd = readtable(cfdfilepath);
    cfd_y = cfd{:,3}.*1000;
    cfd_velo = cfd{:,4};

    %plotting v_bar,vrms, and cfd v_bar for journal paper
    figure;
    t = tiledlayout(1,1);
    ax1 = axes(t);
    plot(ax1,adj.velocity_mean_smooth,adjheights,'b','Linewidth',2.5);
    hold on;
    plot(time_av.cor_velo,time_av.height,'g','Linewidth',2);
    hold on;
    plot(ax1,adj.velocity_mean_smooth-adj.mean_uncertainty_combined,adjheights,':b','Linewidth',1.5);
    hold on;
    plot(ax1,cfd_velo,cfd_y,'k','Linewidth',2);
    hold on;
    plot(ax1,adj.velocity_mean_smooth+adj.mean_uncertainty_combined,adjheights,':b','Linewidth',1.5);
    ax1.XColor = 'b';
    ax1.YColor = 'k';
    ylim(ax1,[0, 10])
    xlim(ax1,[0, 1000])
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times
    xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
    ylabel('Height above the surface [mm]');
    grid on;

    ax2 = axes(t);
    plot(ax2,adj.velocity_rms_smooth,adjheights,'r','Linewidth',2.5);
    hold on;
    plot(ax2,adj.velocity_rms_smooth+adj.rms_uncertainty_combined_perc,adjheights,':r','Linewidth',1.5);
    hold on;
    plot(ax2,adj.velocity_rms_smooth-adj.rms_uncertainty_combined_perc,adjheights,':r','Linewidth',1.5);
    ax2.XColor = 'r';
    ax2.YColor = 'r';
    ax2.XAxisLocation = 'top';
    ax2.YAxis.Visible = 'off';
    ax2.Color = 'none';
    ax1.Box = 'off';
    ax2.Box = 'off';
    ylim(ax2,[0, 10])
    xlim(ax2,[0, 100])
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times
    xlabel('$V_{RMS}$ [m/s]','Interpreter','Latex')
    grid on;


    figure;
    plot(adj.velocity_mean_smooth,adjheights,'b','Linewidth',2.5);
    hold on;
    plot(adj.velocity_mean_smooth-adj.mean_uncertainty_combined,adjheights,':b','Linewidth',1.5);
    hold on;
    plot(adj.velocity_rms_smooth,adjheights,'r','Linewidth',2.5);
    hold on;
    plot(adj.velocity_rms_smooth+adj.rms_uncertainty_combined,adjheights,':r','Linewidth',1.5);
    hold on;

    plot(cfd_velo,cfd_y,'k','Linewidth',2);
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times

%     legend('Mean Velocity','95% CI','RMS Velocity','95% CI','RANS CFD Mean Velocity');
    legend('$\bar{V}$','$\bar{V}$ CI','$V_{RMS}$','$V_{RMS}$ CI','$\bar{V}$ CFD','Interpreter','latex') % R2018b and later
    hold on;
    plot(adj.velocity_mean_smooth+adj.mean_uncertainty_combined,adjheights,':b','Linewidth',1.5);
    hold on;
    plot(adj.velocity_rms_smooth-adj.rms_uncertainty_combined,adjheights,':r','Linewidth',1.5);

%freestream comparison
FLEET_binary = (adjheights>8)&(adjheights<10);
FLEET_freestream = mean(adj.velocity_mean_smooth(FLEET_binary));
CFD_freestream = mean(cfd_velo(150:end));
cfdfleet_relerr = abs(FLEET_freestream-CFD_freestream)./CFD_freestream.*100;

if run_i >1 %not laminar
    %% Inner variable plotting
    cfd_trim_ub = 201;  %154
    cfd_y_plus = cfd{1:cfd_trim_ub,1};
    cfd_u_plus = cfd{1:cfd_trim_ub,2};
    cfd_uvd_plus = cfd{1:cfd_trim_ub,5};
    
    %Calculating inner variables
    Mach = 5.73;
    T_inf = 56; %K, static temp
    P_inf = 370; %Pa, static pressure
    U_inf = 875; %m/s
    y = adjheights;
    u = adj.velocity_mean_smooth;
    [u_plus,u_vd_plus,y_plus,u_plus_c,u_vd_plus_c,y_plus_c] = InnerVariableCalculator(u,y,Mach,T_inf,P_inf,U_inf);
    cfd_y = flipud(cfd{:,3}.*1000);
    cfd_velo = flipud(cfd{:,4});
    U_inf = max(cfd_velo);
    [u_plus_cfd,u_vd_plus_cfd,y_plus_cfd,u_plus_c_cfd,u_vd_plus_c_cfd,y_plus_c_cfd] = InnerVariableCalculator(cfd_velo,cfd_y,Mach,T_inf,P_inf,U_inf);


    %Ideal inner variable scaling
    y_ideal = logspace(1,2);
    xi = 0.41;
    C = 5;
    u_ideal = (1/xi).*log(y_ideal)+C;
        %viscous sublayer
    y_visc = logspace(-1,1);
    u_visc = y_visc;
    
    %plotting innver variable cfd
    figure;
    semilogx(cfd_y_plus,cfd_u_plus,'r','Linewidth',3);
    hold on;
    semilogx(y_plus,u_plus,'b','Linewidth',3);
%     hold on;
%     semilogx(y_plus,u_vd_plus,'b','Linewidth',3);
    hold on;
    semilogx(y_plus_c,u_plus_c,':b','Linewidth',3);
    hold on;
    semilogx(y_ideal,u_ideal,'k','Linewidth',2);
    hold on;
    semilogx(y_visc,u_visc,'--k','Linewidth',2);
    hold on;
    semilogx(y_plus_c_cfd,u_plus_c_cfd,'g','Linewidth',2);

    legend('CFD','FLEET - Fitted Wall Shear','FLEET - Clauser Chart','Logarithmic Law of the Wall (0.41,5)','Viscous Sublayer','dup CFD')
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times
    xlim([0.4,205])
    ylim([0,20])
    grid on;
    xlabel('y^+')
    ylabel('u^+_e_f_f')
    % title('FLEET vs CFD Law of the Wall Profiles')
end

%% Replacing mean velocity below some 3mm height with the time-averaged fit
mod_velo = adj.velocity_mean_smooth(adj_binary);
mod_height = adj.velocimetry_geometricloc(adj_binary,5);
replace_heights_time_ac = mod_height<3;
replace_heights_time_av = time_av.height<3;
mod_velo(replace_heights_time_ac) = time_av.cor_velo(replace_heights_time_av);

%% Output table
time_av_meanvelo = time_av.cor_velo;
Mean_Velocity = adj.velocity_mean_smooth(adj_binary);
Mean_Velocity_u = adj.mean_uncertainty_combined(adj_binary);
RMS_Velocity = adj.velocity_rms_smooth(adj_binary);
RMS_Velocity_u = adj.rms_uncertainty_combined(adj_binary);
Signal_to_Noise_Ratio = adj.mean_SNR(adj_binary);

streamwise = adj.velocimetry_geometricloc(adj_binary,1);
streamwise_u = adj.velocimetry_geometricloc(adj_binary,2);
spanwise = adj.velocimetry_geometricloc(adj_binary,3);
spanwise_u = adj.velocimetry_geometricloc(adj_binary,4);
height = adj.velocimetry_geometricloc(adj_binary,5);
height_u = adj.velocimetry_geometricloc(adj_binary,6)+2.*( adj.velocimetry_geometricloc(10,5)-adj.velocimetry_geometricloc(11,5));

FLEETResultsTable = table(Mean_Velocity,Mean_Velocity_u,RMS_Velocity,RMS_Velocity_u,Signal_to_Noise_Ratio,...
            streamwise,streamwise_u,spanwise,spanwise_u,height,height_u,mod_velo);

%% Save the adjusted data to another file
output_filename = strcat("CorrectedData/CorrectedData_Run",num2str(run_i),"_Images",num2str(maxrealimages));
save(output_filename,'adj','FLEETResultsTable')

% %% FLEET vs DAQ Plotting
% %fleet
% time_accurate_freestream_velo;
% fleetoffset = 8;
% timestep = 0.001;
% timeaverage = 100;
% fleet_time = fleetoffset:timestep:((length(time_accurate_freestream_velo).*timestep)+fleetoffset-timestep);
% time_accurate_freestream_velo(isnan(time_accurate_freestream_velo)) = mean(time_accurate_freestream_velo(~isnan(time_accurate_freestream_velo)));
% time_accurate_freestream_velo_float = movmean(time_accurate_freestream_velo,timeaverage);
% time_accurate_freestream_velo_unc = mean(adj.mean_uncertainty_heightinc(frs_binary)).*ones(size(time_accurate_freestream_velo_float));
% 
% %DAQ
% ACE_timedep_names = "C:\Users\clark\Documents\GitHub\BL_FLEET\ACEDAQPlotter\Run 5246 (ACE).xlsx";
% DAQ_time = readmatrix(ACE_timedep_names,'Sheet','Reduced Data','Range','A2:A1000'); %seconds
% DAQ_velo = readmatrix(ACE_timedep_names,'Sheet','Reduced Data','Range','AW2:AW1000');
% DAQ_time = (0.1:0.1:(length(DAQ_velo))*(1/10))';
% 
% figure;
% plot(DAQ_time,DAQ_velo,'r','Linewidth',3);
% hold on;
% plot(fleet_time,time_accurate_freestream_velo_float,'k','Linewidth',3);
% hold on;
% plot(fleet_time,time_accurate_freestream_velo_float-time_accurate_freestream_velo_unc,'k','Linewidth',1);
% hold on;
% plot(fleet_time,time_accurate_freestream_velo_float+time_accurate_freestream_velo_unc,'k','Linewidth',1);
% 
% xlim([5,20]);
% ylim([825,950]);
% grid on;    
% xlabel('Time [s]')
% ylabel('$\bar{V}_\infty$ [m/s]','Interpreter','Latex')
% legend('ACE DAQ Freestram Velocity','FLEET Freestram Velocity','FLEET 95% CI');
% set(gca,'FontSize', 20);
% set(gca,'fontname','times')  % Set it to times
% 
% %find mean for direct comparison
% fleet_time_avg = mean(time_accurate_freestream_velo);
% DAQ_time_binary = (DAQ_time>min(fleet_time))&(DAQ_time<max(fleet_time));
% DAQ_velo_avg = mean(DAQ_velo(DAQ_time_binary));
% 
% relerr = abs(DAQ_velo_avg-fleet_time_avg)./DAQ_velo_avg;
% disp(['Relative Error in Freestream Velocimetry is ',num2str(relerr.*100),'%'])
% relunc = mean(time_accurate_freestream_velo_unc)./DAQ_velo_avg;
% disp(['Relative Uncertainty in Freestream Velocimetry is ',num2str(relunc.*100),'%'])
% 
