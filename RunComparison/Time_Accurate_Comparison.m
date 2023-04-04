clear all;close all;clc;

%% Time Average Comparison
%This script plots the time-averaged velocity and associated uncertainty
%from the fit to the time-averaged images. 

%% Overarching Variables
folder_path                 = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\ProcessedData\"; 
file_partial_name           = "FullData_Run";
file_partial_name_synth     = "FullData_Synth_Run";
Run_Conditions_filepath     = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/BLFLEETRunConditions.mat';   %stuff like gates and delays
Resolution_filepath         = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/RefData.mat';                    %resolution
ACE_On_Condition_filepath   = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/ACE_Data.mat';                    %resolution
near_wall_folder_path       = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\RunComparison\Near_Wall_Uncertainty_Data\";
near_wall_file_name         = "Time_Average_Fit_RR";
near_wall_folder_path_5479  = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\RunComparison\Near_Wall_Uncertainty_Data_5479\";
inlet_cfd_filepath          = "Sept2022_RANS_CFD.xlsx";
DNS_filepath                = "CFDDataForPlot.xlsx";

    load(Run_Conditions_filepath);
    load(Resolution_filepath);
    load(ACE_On_Condition_filepath);

runs_list = 1:14;
lens_standoff_dist = 300; %mm
height_focus = 3; %focusing around 3mm from the surface
decay_folderpath = "C:\Users\clark\Documents\GitHub\BL_FLEET_2.0\RawDataProcessing\ProcessedData\";
decay_filepaths = [ "ProcessedData_Synth_Run1_Images1500.mat";...
                    "ProcessedData_Synth_Run2_Images300.mat";...
                    "ProcessedData_Synth_Run3_Images300.mat"];
FLEET_July_2022_filepath = "C:\Users\clark\Documents\Grad_Research\Data\FLEET_July2022\Data_Compiled.xlsx";
dist_LE_Trips = 54.065; %mm, distance from the leading edge to the center of the tripping element

colorlist = ['r','b','m','g'];
cfd_colorlist = ['k',':k'];
PB_repeat   = [5484,5485,5486];
SRA_repeat  = [5479,5483,5495];
PB_BL       = [5491,5490,PB_repeat];
PB_BL_sp    = [5491,5490,5486];
SRA_BL      = [5492,5493,5494,SRA_repeat];
SRA_BL_sp   = [5493,5494,5495];
PB_downstm  = [5487,5488,5489]; %C, R, C-downstream
all_runs    = [5479,5483:5495];
locations_list = [3;3;3;3;3;4;5;6;2;1;1;1;2;3];

%% Getting more complex input variables
    %load in resolution and gate/delay pairs
        load(Run_Conditions_filepath);
        load(Resolution_filepath);
    %find filepaths
            a=dir(strcat(folder_path,file_partial_name, "*.mat"));
            all_filenames = strings(length(uniqueruns),1);
            for i = 1:length(a)
                all_filenames(i) = a(i).name;
            end
            %get only the appropriate runs
            ordered_filepaths = strings(length(uniqueruns),1);
            for i = 1:length(ordered_filepaths)
                filename_exp = strcat(file_partial_name,num2str(runs_list(i)),"_");
                which_run_binary = contains(all_filenames,filename_exp);
                ordered_filepaths(i) = all_filenames(which_run_binary);
            end

    %Load in variables
    c_centroids         = cell(length(runs_list),1);
    c_centroid_erros    = cell(length(runs_list),1);
    c_heights           = cell(length(runs_list),1);
    c_SNRs              = cell(length(runs_list),1);
    c_bounds            = cell(length(runs_list),2);
    c_gates             = cell(length(runs_list),1);
    c_delays            = cell(length(runs_list),1);
    c_resolutions       = cell(length(runs_list),1);
    c_downstream_locs   = cell(length(runs_list),1);
    c_downstream_locs_unc= cell(length(runs_list),1);
    c_spanwise_locs     = cell(length(runs_list),1);
    c_spanwise_locs_unc = cell(length(runs_list),1);
    c_wall_location_unc = cell(length(runs_list),1);
    c_decay             = cell(length(runs_list),1);
    c_num_ims           = cell(length(runs_list),1);

    %calculated/derived variables
    c_velocities                = cell(length(runs_list),1);
    c_rms                       = cell(length(runs_list),1);
    c_heights_adjusted          = cell(length(runs_list),1);
    c_velo_bounds               = cell(length(runs_list),1);  
    c_SNR_binary_filt           = cell(length(runs_list),1);

    %Uncertainties
    unc_velo_centroid_fit       = cell(length(runs_list),1);    %from fitting
    unc_velo_resolution         = cell(length(runs_list),1);    %from fitting
    unc_rms_centroid_fit        = cell(length(runs_list),1);    %from fitting
    unc_rms_resolution          = cell(length(runs_list),1);    %from fitting
    unc_height_resolution       = cell(length(runs_list),1);   
    unc_height_wall_location    = cell(length(runs_list),1); 
    unc_velo_emission_decay     = cell(length(runs_list),1); 
    unc_rms_emission_decay      = cell(length(runs_list),1); 
    unc_velo_magnification      = cell(length(runs_list),1); 
    unc_rms_magnification       = cell(length(runs_list),1); 
    unc_height_magnification    = cell(length(runs_list),1); 
    unc_height_inclination      =cell(length(runs_list),1); 
    unc_velo_span_point         = cell(length(runs_list),1);
    unc_rms_span_point          = cell(length(runs_list),1);
    unc_velo_nearwall           = cell(length(runs_list),1);
    unc_velo_RANDOM             = cell(length(runs_list),1);
    unc_velo_SYSTEMATIC         = cell(length(runs_list),1);
    unc_velo_COMBINED           = cell(length(runs_list),1);
    unc_rms_RANDOM             = cell(length(runs_list),1);
    unc_rms_SYSTEMATIC         = cell(length(runs_list),1);
    unc_rms_COMBINED           = cell(length(runs_list),1);
    unc_height_SYSTEMATIC       = cell(length(runs_list),1);
    unc_velo_heights            = cell(length(runs_list),1); 

    %Inner Variables
    inner_velo_fleet    = cell(length(runs_list),1); 
    inner_y_fleet       = cell(length(runs_list),1);
    inner_velo_cfd      = cell(length(runs_list),1);
    inner_y_cfd         = cell(length(runs_list),1); 

    %decay constant data sets
    c_taus          = cell(length(decay_filepaths),1);
    c_tau_errors    = cell(length(decay_filepaths),1);

    %synthetic data variables
    s_velocities             = cell(length(runs_list),1); 
    s_rms                    = cell(length(runs_list),1); 
    s_rms_uncertainty        = cell(length(runs_list),1); 
    s_SNRs                   = cell(length(runs_list),1); 
    s_decay                  = cell(length(runs_list),1); 
    s_tau_input              = cell(length(runs_list),1); 
    s_velo_corr              = cell(length(runs_list),1); 
    s_heights                = cell(length(runs_list),1); 
    unc_rms_synthetic        = cell(length(runs_list),1); 

    %save out data into an excel file
        %headers
        table_headers = {'Mean Velocity [m/s]','Uncertainty in Mean Velocity [m/s]',...
                        'RMS Velocity [m/s]','Uncertainty in RMS Velocity [m/s]','Height above the Test Article [mm]',...
        	            'Uncertainty in Height above the Test Article [mm]','Location Downstream the Leading Edge [mm]',...
                    	'Uncertainty in Downstream Location  [mm]','Spanwise Location [mm]','Uncertainty in Spanwise Location [mm]'};
        excel_filepath = "FLEET_Sept2022.xlsx";

%% Loading in Data from the Sept 2022 Campaign
for i= 1:length(runs_list)
    load(fullfile(folder_path,ordered_filepaths(i)))

    if i ~=14
    c_velocities{i}             = velocity_mean;
    c_rms{i}                    = velocity_rms;
    c_heights{i}                = velocimetry_geometricloc(:,5)';
    c_gates{i}                  = Gates(i,:);
    c_delays{i}                 = Delays(i,:);
    if i~=11
    c_SNRs{i}                   = mean_SNR;
    else
    c_SNRs{i}                   = mean_SNR.^(2/3);
    end
    c_resolutions{i}            = pixel_um_resolution(i,:);
    c_wall_location_unc{i}      = zero_height_ref_unc;
    c_decay{i}                  = tau_fit;
    unc_velo_centroid_fit{i}    = velocity_mean_r;
    unc_velo_resolution{i}      = velocity_mean_s;
    unc_rms_centroid_fit{i}     = velocity_rms_r;
    unc_rms_resolution{i}       = velocity_rms_s;
    c_num_ims{i}                = invar_image_compare;

    else %run 14
        height_binary               = velocimetry_geometricloc(:,5)'<8.5; 
        c_velocities{i}             = velocity_mean(height_binary);
        c_rms{i}                    = velocity_rms(height_binary);
        c_heights{i}                = velocimetry_geometricloc(height_binary',5)';
        c_gates{i}                  = Gates(i,:);
        c_delays{i}                 = Delays(i,:);
        c_SNRs{i}                   = mean_SNR(height_binary);
        c_resolutions{i}            = pixel_um_resolution(i,:);
        c_wall_location_unc{i}      = zero_height_ref_unc;
        c_decay{i}                  = tau_fit(height_binary);
        unc_velo_centroid_fit{i}    = velocity_mean_r(height_binary);
        unc_velo_resolution{i}      = velocity_mean_s(height_binary);
        unc_rms_centroid_fit{i}     = velocity_rms_r(height_binary);
        unc_rms_resolution{i}       = velocity_rms_s(height_binary);
        c_num_ims{i}                = invar_image_compare;
    end

end

%% Load in CFD Data
[aw,no_aw] = LoadCFD_Data_Sept2022(inlet_cfd_filepath);

%% Loading in Data for near-wall uncertainty 
[near_wall_fitobject] = NearWall_Uncertainty_Calculator(near_wall_folder_path,near_wall_file_name,Gates,Delays,pixel_um_resolution);
[near_wall_fitobject_5479] = NearWall_Uncertainty_Calculator(near_wall_folder_path_5479,near_wall_file_name,Gates,Delays,pixel_um_resolution);

%% Load in Synthetic Data
%find filepaths
        a=dir(strcat(folder_path,file_partial_name_synth, "*.mat"));
        all_filenames = strings(length(a),1);
        for i = 1:length(a)
            all_filenames(i) = a(i).name;
        end
        %get only the appropriate runs
        ordered_filepaths = strings(length(uniqueruns),1);
        for i = 1:length(ordered_filepaths)
            filename_exp = strcat(file_partial_name_synth,num2str(runs_list(i)),"_");
            which_run_binary = contains(all_filenames,filename_exp);
            ordered_filepaths(i) = all_filenames(which_run_binary);
        end

for j = 1:length(ordered_filepaths)
     
    load(fullfile(folder_path,ordered_filepaths(j)));

    s_velocities{j}             = velocity_mean;
    s_rms{j}                    = velocity_rms;
    s_SNRs{j}                   = mean_SNR;
    s_rms_uncertainty{j}        = velocity_rms_s;
    s_decay{j}                  = tau_fit;
    s_tau_input{j}              = synth_input_tau_fit;
    s_velo_corr{j}              = (nondim_velo_error).*synth_input_velocity_mean-synth_input_velocity_mean;
    s_heights{j}                = velocimetry_geometricloc(:,5);
    mean_decay(j)               = mean(synth_input_tau_fit);
    synth_mean_velo(j)          = mean(synth_input_velocity_mean);
end

%% Loading in DNS simulation (very different run conditions, just for a rough comparison)
DNS_V = readtable(DNS_filepath,'Range','E2:E202');
DNS_RMS = readtable(DNS_filepath,'Range','G2:G202');
DNS_h = readtable(DNS_filepath,'Range','J2:J202');
DNS_yplus = readtable(DNS_filepath,'Range','B2:B202');
DNS_uplus = readtable(DNS_filepath,'Range','F2:F202');
DNS_V = DNS_V{:,:};
DNS_RMS = DNS_RMS{:,:};
DNS_h = DNS_h{:,:};
DNS_yplus = DNS_yplus{:,:};
DNS_uplus = DNS_uplus{:,:};

%% Imprecision Correction using the RMS correction 
RMS_precision_adj = [10;15;-8;-15;-15;-5;0;0;0;0;0;-15;0;0];
for j= 1:length(runs_list)

    %measured data
    meas_rms =  c_rms{j};
    meas_snr =  c_SNRs{j};
    meas_rms_unc = unc_rms_resolution{j};
    meas_heights = transpose(c_heights{j});

    %fit synthetic data
    use_synth = 14;
    snr = s_SNRs{use_synth};
    rms = s_rms{use_synth};
    unc = s_rms_uncertainty{use_synth};
    heights = s_heights{use_synth};
    heights_binary = and((heights>1),snr>6);%mm
    snr = snr(heights_binary);
    rms = rms(heights_binary)+RMS_precision_adj(j);
    unc = unc(heights_binary);
    synthrmsfit=fit(snr,rms,'power2'); 
    fit_conf = confint(synthrmsfit);
    synth_minfit = @(x) fit_conf(1,1).*x.^(fit_conf(1,2))+fit_conf(1,3);
    synth_maxfit = @(x) fit_conf(2,1).*x.^(fit_conf(2,2))+fit_conf(2,3);
    snr_stepper = linspace(min(snr),max(snr)+100,300);
          
    %calculate the corrected rms velocity
    rms_imprecision = synthrmsfit(meas_snr);
    min_rms_imprecision = synth_minfit(meas_snr);
    max_rms_imprecision = synth_maxfit(meas_snr);
    synthfit_conf = (max_rms_imprecision-min_rms_imprecision)./2;

    for i = 1:length(rms_imprecision)
        if rms_imprecision(i) > meas_rms(i)
        rms_imprecision(i) = meas_rms(i);
        end
    end

    rms_real        = sqrt(meas_rms.^2-rms_imprecision.^2);
    rms_real_unc    =  sqrt((meas_rms_unc.^2)+(mean(unc)).^2-meas_rms_unc.*mean(unc));
    synthfit_conf(isnan(rms_real)) = NaN;
%     %original from before 1/20

    %Only plot heights below a specified height where uncertainty is large
    %(crosses y axis)
    min_possible_rms = rms_real-rms_real_unc;
    less_zero_binary = min_possible_rms>0;
    [r,c]=find(less_zero_binary);
    [r,ir]=unique(r,'first');
    max_height = meas_heights(cumsum(cumsum(less_zero_binary)) == 1);
    height_bounding_binary = meas_heights<=max_height;
    
        %only use heights below max height
        rms_real(~height_bounding_binary) = NaN;

    %specialty trimming because I don't want to write a program to pick out
    %every time the uncertainty bar crosses zero and get fancy
    if j == 1 %5479
     height_bounding_binary = meas_heights<=6.5;
        rms_real(~height_bounding_binary) = NaN;
    elseif j== 10 %5491
    height_bounding_binary = meas_heights<=4.56;
        rms_real(~height_bounding_binary) = NaN;
    elseif j==11 %5492
    height_bounding_binary = meas_heights>=3.15;
        rms_real(~height_bounding_binary) = NaN;
    elseif j==12 %5493
    height_bounding_binary = meas_heights>=1.02;
        rms_real(~height_bounding_binary) = NaN;
    elseif j== 13 %5494
    height_bounding_binary = meas_heights<=5.5;
        rms_real(~height_bounding_binary) = NaN;
    elseif j== 14 %5494
    height_bounding_binary = meas_heights<=6.5;
        rms_real(~height_bounding_binary) = NaN;
    end

    %aSaving the correction
    c_rms{j} = rms_real;
    unc_rms_synthetic{j} = synthfit_conf;

if j==9

    figure;
    subplot(1,2,1);
    scatter(snr,rms,'r');
    hold on;
    plot(snr_stepper,synthrmsfit(snr_stepper),'r','Linewidth',2);
%     errorbar(snr,rms,unc,'b');
    scatter(meas_snr,meas_rms,'b','filled');
    scatter(meas_snr,rms_real,'m');
%     errorbar(meas_snr,rms_real,rms_real_unc,'k');
    xlabel('SNR [-]');
    ylabel('$V_{RMS} [m/s]$','Interpreter','Latex')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    grid on;
%     title('Imprecision estimate')
    ylim([20,max(rms)+30])
    legend('Synthetic','Fit to Synthetic','Measured','Corrected RMS Velocity');
    xlim([0,max(meas_snr)+10])
    set(gca,'FontSize', 18);
    set(gca,'fontname','times')  % Set it to times

    subplot(1,2,2);
    plot(rms_real,meas_heights,'r','Linewidth',2);
    xlabel('$Corrected V_{RMS} [m/s]$','Interpreter','Latex')
    ylabel('Height above surface [mm]');
    title(['Run ',num2str(uniqueruns(j))]);
    xlim([0,max(rms_real)+25])

    figure;
    scatter(snr,rms,'r');
    hold on;
    plot(snr_stepper,synthrmsfit(snr_stepper),'r','Linewidth',2);
%     errorbar(snr,rms,unc,'b');
    scatter(meas_snr,meas_rms,'b','filled');
    scatter(meas_snr,rms_real,'m');
%     errorbar(meas_snr,rms_real,rms_real_unc,'k');
    xlabel('SNR [-]');
    ylabel('$V_{RMS} [m/s]$','Interpreter','Latex')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    grid on;
%     title('Imprecision estimate')
    ylim([20,max(rms)+30])
%     legend(["Synthetic: $V_{RMS,\sigma-err}(SNR)$","Fit: $V_{RMS,\sigma-err}(SNR)$",...
%         "Measured: $V_{RMS}$","Real Fluctuations: $V_{RMS,\sigma-corr}$"],'Interpreter','latex');
      legend(["Synthetic: $V_{RMS,err}(SNR)$","Fit: $V_{RMS,err}(SNR)$",...
        "Measured: $V_{RMS,m}$","Real Fluctuations: $V_{RMS}$"],'Interpreter','latex');
    xlim([0,max(meas_snr)+10])
    set(gca,'FontSize', 18);
    set(gca,'fontname','times')  % Set it to times

    figure;
    plot(meas_snr,meas_heights,'k','Linewidth',2);
    xlabel('SNR [-]');
    ylabel('Height above surface [mm]')
    set(gca,'FontSize', 18);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,max(meas_snr)+10])
    ylim([0,12]);
    grid on;


figure;
plot(meas_rms,meas_heights,'b','Linewidth',2);
hold on;
plot(synthrmsfit(meas_snr),meas_heights,':r','Linewidth',2);
plot(rms_real,meas_heights,'m','Linewidth',2);

xlabel('$V_{RMS}$ [m/s]','Interpreter','latex');
ylabel('Height above surface [mm]')
ylim([0,10]);
set(gca,'FontSize', 18);
set(gca,'fontname','times')  % Set it to times
legend(["Measured: $V_{RMS,m}$","Imprecision: $V_{RMS,err}(SNR)$","Real Fluctuations: $V_{RMS}$"],'Interpreter','latex');
grid on;
ylim([0,12]);

% 
% figure;
% plot(meas_snr,meas_heights,'b','Linewidth',2);
% xlabel('SNR');
% ylabel('Height above surface [mm]')
% ylim([0,10]);
% xlim([0,max(meas_snr)+10])
% set(gca,'FontSize', 18);
% set(gca,'fontname','times')  % Set it to times
% ylim([0,12]);



disp('wait');
close all;
end

end

%% Correcting for Emission Decay
decay_offset = cell(length(runs_list),1); 
decay_offset_freestream = cell(length(runs_list),1);
%get synthetic information from previous simulations
for i = 1:length(runs_list)
    %Synthetic stuff
            %fit the decay constant vs error of the synthetic data
            %load synthetic data
            tau_fit             = s_decay{i};
            synth_input_tau_fit = s_tau_input{i};
            velo_error          = s_velo_corr{i};

            %fit the synthetic data
            s_tau_lin = polyfit(synth_input_tau_fit,velo_error,1);

            %plot the quality of the fit to the synthetic data
            plot_tau = linspace(min(synth_input_tau_fit),max(synth_input_tau_fit),100);
            velo_fit_synth = polyval(s_tau_lin,plot_tau);

%                 figure;
%                 hold on;
%                 scatter(synth_input_tau_fit,velo_error,'r');
%                 hold on;
%                 plot(plot_tau,velo_fit_synth,'k','Linewidth',2);
%                 xlabel('Decay Constant [ns]');
%                 ylabel('Velocity');
%                 title('Fitting decay constant vs velocity error in the syntehtic dataset')
            
        %Finding a linear-approximation of the decay constant in the real
        %data
                % Loading in real data
                    tau_run = c_decay{i};
                    heights_run = c_heights{i};
                    SNR_run = c_SNRs{i}';
        
                %use the snr binary filter to only get data with a good SNR above
                %the cutoff
                    switchover_height_mm = 2;  %mm
                    tau_filt = heights_run>=switchover_height_mm;
                    
                %filter
                    tau_run_f = tau_run(tau_filt);
                    heigths_run_f = heights_run(tau_filt);
        
                %fit
                    tau_fit_run = polyfit(heigths_run_f,tau_run_f,1);
                    lin_fit_taus = polyval(tau_fit_run,heights_run);       

                %plotting the linear fit of the decay constant (no velocity
                %yet)
%                     figure;
%                     hold on;
%                     scatter(tau_run,heights_run,'r');
%                     hold on;
%                     scatter(tau_run_f,heigths_run_f,'b');
%                     plot(lin_fit_taus,heights_run,'k','Linewidth',2);
%                     xlabel('Decay Constant [ns]');
%                     ylabel('Height Above the Surface');
%                     title(['Fitting tau as a line for Run ',num2str(uniqueruns(i))])
% 

        %Using the linear-fit-measured-tau to extrapolate the error in
        %velocity 
        correcting_for_tau_constant = 0;
        decay_offset_freestream{i} = polyval(s_tau_lin,lin_fit_taus-correcting_for_tau_constant);
        decay_offset{i}= (0.75).*(c_velocities{i}./(max(c_velocities{i}))).*decay_offset_freestream{i}';
        tau_unc = [transpose(lin_fit_taus-40),transpose(lin_fit_taus+40)];
        decay_offset_lims =  polyval(s_tau_lin,tau_unc);
        decay_offset_uncertainty = abs(decay_offset_lims(:,1)-decay_offset_lims(:,2))/2;

        disp(['Decay correction for Run ',num2str(uniqueruns(i)),' is ',num2str(mean(decay_offset{i})),...
            ' with bounds of ',num2str(min(decay_offset{i})),' to ', num2str(max(decay_offset{i})),''])

             c_velocities{i} = c_velocities{i} + (decay_offset{i}).*(c_velocities{i}./max(c_velocities{i}));
             unc_velo_emission_decay{i} = decay_offset_uncertainty;
             unc_rms_emission_decay{i} = 0;
% 
%                 figure;
%                 subplot(1,2,1);
%                 scatter(tau_run_f,heigths_run_f,'b');
%                 plot(lin_fit_taus,heights_run,'k','Linewidth',2);
%                 xlabel('Decay Constant [ns]');
%                 ylabel('Height Above the Surface');
%                 title(['Real Tau'])
%                 grid on;
% 
%                 subplot(1,2,2);
%                 plot(decay_offset,heights_run,'k','Linewidth',2);
%                 xlabel('Error in mean velocity');
%                 ylabel('Height Above the Surface');
%                 title(['Real Tau'])
%                 grid on;
end

%% Uncertainty Propogation (wall location, magnification, inclination angle, span angle, near-wall fitting effects)
for i= 1:length(runs_list)

    %     SpanwiseCameraAngle
    %     SpanwiseCameraAngle_unc %second column of SpanwiseCameraAngle

    %wall_loc_unc-------------------------------------------------------------
            heights = c_heights{i}';
            wall_loc_unc = sqrt(c_wall_location_unc{i}.^2+1.^2);
            res = c_resolutions{i};
            [resolution_height_uncertainty,wall_location_height_uncertainty] ...
                = Wall_Resolution_height_Uncertainty_Propogator(res,wall_loc_unc,heights);

            unc_height_resolution{i} = resolution_height_uncertainty;
            unc_height_wall_location{i} = wall_location_height_uncertainty.*ones(size(resolution_height_uncertainty));

 %Inclination angle
        c_heights{i} = c_heights{i}./cosd(inclination_angle(i));
        unc_height_inclination{i} = c_heights{i}'.*(sind(inclination_angle(i))./(cosd(inclination_angle(i).^2))).*deg2rad(inclination_angle_unc(i));

    %magnification and beam angle
        L_hvec = heights-3; 
        angle_beam = 0; %likely beam angle, degrees
        angle_beam_unc = 4; %likely beam angle, degrees

        %thesis definitions
        delta_M = (((L_hvec.*cosd(angle_beam))/(lens_standoff_dist))).*((deg2rad(angle_beam_unc)));
        unc_velo_magnification{i} = c_velocities{i}.*delta_M;
        unc_rms_magnification{i} = c_rms{i}.*delta_M;
        unc_height_magnification{i} = c_heights{i}'.*delta_M;

    %Span angle
        %change systematic errors
            span_angle =  SpanwiseCameraAngle(i,1); %degree
            span_angle_uncertainty = 2*SpanwiseCameraAngle(i,2); %degree
        %correct mean and RMS velocity
            c_velocities{i} = c_velocities{i}./cosd(span_angle);
            c_rms{i} = c_rms{i}./cosd(span_angle);
        %account for uncertainty propogation
            unc_velo_span_point{i} = c_velocities{i}.*(sind(span_angle)./(cosd(span_angle.^2))).*deg2rad(span_angle_uncertainty);
            unc_rms_span_point{i}  = c_rms{i}.*(sind(span_angle)./(cosd(span_angle.^2))).*deg2rad(span_angle_uncertainty);
    %Near-wall effects
    if i==1 %first run, different config
        unc_velo_nearwall{i} = feval(near_wall_fitobject_5479,c_heights{i});
    else
        unc_velo_nearwall{i} = feval(near_wall_fitobject,c_heights{i});
    end

end

%% Combine random and systematic uncertainties
for i = 1:length(runs_list)

    %mean velocity
    unc_velo_RANDOM{i} = sqrt(unc_velo_centroid_fit{i}.^2);
    unc_velo_SYSTEMATIC{i} = sqrt(unc_velo_resolution{i}.^2+unc_velo_emission_decay{i}.^2+unc_velo_magnification{i}.^2 ...
                                +unc_velo_span_point{i}.^2+unc_velo_nearwall{i}.^2);
    unc_velo_COMBINED{i} = sqrt(unc_velo_RANDOM{i}.^2+unc_velo_SYSTEMATIC{i}.^2);

    %rms velocity
    unc_rms_RANDOM{i} = sqrt(unc_rms_centroid_fit{i}.^2);
%     unc_rms_SYSTEMATIC{i} = sqrt(unc_rms_resolution{i}.^2+unc_rms_magnification{i}.^2 ...
%                                 +unc_rms_span_point{i}.^2+unc_rms_synthetic{i}.^2);
    unc_rms_SYSTEMATIC{i} = sqrt(unc_rms_magnification{i}.^2 ...
                                +unc_rms_span_point{i}.^2+unc_rms_synthetic{i}.^2);
    unc_rms_COMBINED{i} = sqrt(unc_rms_RANDOM{i}.^2+unc_rms_SYSTEMATIC{i}.^2);

    %heights
    unc_height_SYSTEMATIC{i} = sqrt(unc_height_resolution{i}.^2+unc_height_wall_location{i}.^2+unc_height_inclination{i}.^2+unc_height_magnification{i}.^2);
    c_heights{i} = c_heights{i}';
end

for i= 1:length(unc_rms_COMBINED)
    allvals = unc_rms_COMBINED{i};
test(i) = mean(allvals(~isnan(allvals)));
end


%% Make a table showing the reletive intensities of each form of uncertainty for height, RMS, and mean velocity
    for m = 1:1
        
        %height
        table_headers_height = {'Resolution','Wall Location','Inclination','BeamPointing'};
        
            for i = 1:length(runs_list)
                h_res(i) = 100.*mean(unc_height_resolution{i})./max(c_heights{i});
                h_wall(i) = 100.*mean(unc_height_wall_location{i})./max(c_heights{i});
                h_inc(i) = 100.*mean(unc_height_inclination{i})./max(c_heights{i});
                h_mag(i) = 100.*mean(unc_height_magnification{i})./max(c_heights{i});
            end
                h_res = mean(h_res);
                h_wall = mean(h_wall);
                h_inc = mean(h_inc);
                h_mag = mean(h_mag);
            
                T_height= table(h_res',h_wall',h_inc',h_mag','VariableNames',table_headers_height);
        
        %mean
        
        perc_max = 90;
        perc_min = 10;
        table_headers_velo = {'Centroids','Resolution','Emission Decay','BeamPointing','SpanwisePointing','Nearwall'};
        
            for i = 1:length(runs_list)
                V_cent(i) = 100.*prctile(unc_velo_centroid_fit{i},perc_max)./max(c_velocities{i});
                V_res(i) = 100.*prctile(unc_velo_resolution{i},perc_max)./max(c_velocities{i});
                V_decay(i) = 100.*prctile(unc_velo_emission_decay{i},perc_max)./max(c_velocities{i});
                V_mag(i) = 100.*prctile(unc_velo_magnification{i},perc_max)./max(c_velocities{i});
                V_span(i) = 100.*prctile(unc_velo_span_point{i},perc_max)./max(c_velocities{i});
                V_wall(i) = 100.*prctile(unc_velo_nearwall{i},98)./max(c_velocities{i});
        %     
        %         V_cent_min(i) = prctile(unc_velo_centroid_fit{i},perc_min)./max(c_velocities{i});
        %         V_res_min(i) = prctile(unc_velo_resolution{i},perc_min)./max(c_velocities{i});
        %         V_decay_min(i) = prctile(unc_velo_emission_decay{i},perc_min)./max(c_velocities{i});
        %         V_mag_min(i) = prctile(unc_velo_magnification{i},perc_min)./max(c_velocities{i});
        %         V_span_min(i) = prctile(unc_velo_span_point{i},perc_min)./max(c_velocities{i});
        %         V_wall_min(i) = prctile(unc_velo_nearwall{i},perc_min)./max(c_velocities{i});
        
        
            end
            
                V_cent = mean(V_cent);
                V_res = mean(V_res);
                V_decay = mean(V_decay);
                V_mag = mean(V_mag);
                V_span = mean(V_span);
                V_wall = mean(V_wall);
        %         V_cent_min = mean(V_cent_min);
        %         V_res_min = mean(V_res_min);
        %         V_decay_min = mean(V_decay_min);
        %         V_mag_min = mean(V_mag_min);
        %         V_span_min = mean(V_span_min);
        %         V_wall_min = mean(V_wall_min);
        %         V_cent = [V_cent,V_cent_min];
        %         V_res = [V_res,V_res_min];
        %         V_decay = [V_decay,V_decay_min];
        %         V_mag = [V_mag,V_mag_min];
        %         V_span = [V_span,V_span_min];
        %         V_wall = [V_wall,V_wall_min];
        
        
            
                T_velo= table(V_cent',V_res',V_decay',V_mag',V_span',V_wall','VariableNames',table_headers_velo);
        
        %RMS
        
        unc_rms_RANDOM{i} = (1/sqrt(c_num_ims{i})).*sqrt(unc_rms_centroid_fit{i}.^2);
            unc_rms_SYSTEMATIC{i} = sqrt(unc_rms_resolution{i}.^2+unc_rms_magnification{i}.^2 ...
                                        +unc_rms_span_point{i}.^2+unc_rms_synthetic{i}.^2);
        
        
        perc_max = 90;
        perc_min = 10;
        table_headers_RMS = {'Centroids','Resolution','Imprecision','BeamPointing','SpanwisePointing'};
        
            for i = 1:length(runs_list)
                norm_rms = 0.12.*max(c_velocities{i});
                VRMS_cent(i) = 100.*prctile(unc_rms_centroid_fit{i},perc_max)./norm_rms;
                VRMS_res(i) = 100.*prctile(unc_rms_resolution{i},perc_max)./norm_rms;
                VRMS_imp(i) = 100.*prctile(unc_rms_synthetic{i},perc_max)./norm_rms;
                VRMS_mag(i) = 100.*prctile(unc_rms_magnification{i},perc_max)./norm_rms;
                VRMS_span(i) = 100.*prctile(unc_velo_span_point{i},perc_max)./norm_rms;
        
        
            end
%             
                VRMS_cent = mean(VRMS_cent);
                VRMS_res = mean(VRMS_res);
                VRMS_imp = mean(VRMS_imp);
                VRMS_mag = mean(VRMS_mag);
                VRMS_span = mean(VRMS_span);
%             
                T_RMS= table(VRMS_cent',VRMS_res',VRMS_imp',VRMS_mag',VRMS_span','VariableNames',table_headers_RMS);

    end

%% Smoothing the mean velocity and uncertainty for plotting
for i = 1:length(runs_list)
    r_velo      = c_velocities{i};
    r_velo_err  = unc_velo_COMBINED{i};
    r_height    = c_heights{i};
    r_rms      = c_rms{i};
    r_rms_err  = unc_rms_COMBINED{i};
    notnan_rms = ~isnan(r_rms);

        res = c_resolutions{i};
        pix_wid = 0.75;
        wid = pix_wid.*res(1)./1000; %pixels into mm
    smooth_velo     = smooth(r_height,r_velo,wid,'lowess');
    smooth_velo_err = smooth(r_height,r_velo_err,wid,'lowess');

%     %plotting
%         figure(1);
%         subplot(1,2,1);
%         plot(c_velocities{i},c_heights{i},'b','Linewidth',2);
%         hold on;
%         plot(c_velocities{i}-r_velo_err,c_heights{i},'k');
%         plot(c_velocities{i}+r_velo_err,c_heights{i},'k');
%     
%         plot(smooth_velo,c_heights{i},'r','Linewidth',2);
%         plot(smooth_velo-smooth_velo_err,c_heights{i},':r');
%         plot(smooth_velo+smooth_velo_err,c_heights{i},':r');
%         title(['Velocity Smoothing for Run', num2str(uniqueruns(i))]);
%         xlabel('Velocity [m/s]');
%         ylabel('Height above the surface [mm]');
%         set(gca,'FontSize', 15);
%         set(gca,'fontname','times')  % Set it to times
%         xlim([0,900]);
%         ylim([0,max(r_height)+1])
%         hold off;
%     
%         subplot(1,2,2);
%         plot(c_velocities{i},c_heights{i},'b','Linewidth',2);
%         xlim([0,1000]);
%         ylim([0,max(c_heights{i})])
%         xlabel('Mean Velocity')
%         ylabel('Height above surface [mm]')
%         grid on;
%         set(gca,'FontSize', 15);
%         set(gca,'fontname','times')  % Set it to times

    %replace only data above where I expect wall effects
    wall_effect_height = 0.5; %mm
    replace_binary = r_height>wall_effect_height;
    r_velo(replace_binary)      = smooth_velo(replace_binary);
    r_velo_err(replace_binary)  = smooth_velo_err(replace_binary);

    c_velocities{i}             = r_velo;
    unc_velo_COMBINED{i}        = r_velo_err;
end

%% Propogating uncertainties in height into uncertainty in the mean velocity (mostly for plotting purposes)
unc_velo_COMBINED_table = unc_velo_COMBINED;
unc_rms_COMBINED_table = unc_rms_COMBINED;

for i = 1:length(runs_list)
velo_unc        =  unc_velo_COMBINED{i};
rms_unc         =  unc_rms_COMBINED{i};
height_unc      =  unc_height_SYSTEMATIC{i}; 
velo            =  c_velocities{i};
rms             =  c_rms{i};
heights         =  c_heights{i};

        %move height up and down
            heights_dn = heights-height_unc;
            heights_up = heights+height_unc;
        %interpolate results if moved up or down
            velo_dn = interp1(heights_dn,velo,heights,'spline');
            velo_up = interp1(heights_up,velo,heights,'spline');
            rms_dn = interp1(heights_dn,rms,heights,'spline');
            rms_up = interp1(heights_up,rms,heights,'spline');
        %provide result
            uncertainty_mean_velo_wall_loc = (abs(velo_dn-velo)+abs(velo_up-velo))/2;
            uncertainty_rms_wall_loc       = (abs(rms_dn-velo)+abs(rms_up-velo))/2;
            unc_velo_heights{i}     = uncertainty_mean_velo_wall_loc;
            unc_velo_COMBINED{i}    = sqrt(unc_velo_COMBINED{i}.^2+unc_velo_heights{i}.^2);
            unc_rms_COMBINED{i}     = sqrt(unc_rms_COMBINED{i}.^2+unc_velo_heights{i}.^2);


            %plot
%                 figure;
%                 subplot(1,2,1);
%                 plot(velo,heights,'b','Linewidth',2);
%                 hold on;
%                 plot(velo,heights_dn);
%                 plot(velo,heights_up);
%                 plot(velo_dn,heights,':r','Linewidth',1);
%                 plot(velo_up,heights,':r','Linewidth',1);
%                     plot(velo-uncertainty_mean_velo_wall_loc,heights,':b','Linewidth',2);
%                     plot(velo+uncertainty_mean_velo_wall_loc,heights,':b','Linewidth',2);
%                 grid on;
%                 xlabel('Mean Velocity [m/s]')
%                 ylabel('Height above the surface [mm]')
%                 title('Uncertainty due to wall location Finding')
% 
%                 subplot(1,2,2);
%                 plot(uncertainty_mean_velo_wall_loc,heights,'b','Linewidth',2);
%                 xlim([0,50]);
%                 xlabel('Velocity Uncertainty [m/s]')
%                 title('Uncertainty due to wall location Finding')

end

%% Saving data into a table for CFD Validation
for i = 1:length(runs_list)

    velo            = round(c_velocities{i},1);
    velo_unc        = round(unc_velo_COMBINED_table{i},1);
    rms             = round(c_rms{i},1);
    rms_unc         = round(unc_rms_COMBINED_table{i},2);
    heights         = round(c_heights{i},3);
    heights_unc     = round(unc_height_SYSTEMATIC{i},3);
    down_loc        = downstream_loc(i).*ones(size(velo));
    down_loc_unc    = downstream_loc_unc(i).*ones(size(velo));
    span_loc        = spanwise_loc(i).*ones(size(velo));
    span_loc_unc    = spanwise_loc_unc(i).*ones(size(velo));

    T = table(velo,velo_unc,rms,rms_unc,heights,heights_unc,down_loc,down_loc_unc,span_loc,span_loc_unc,'VariableNames',table_headers);
    sheet_name = ['Run ',num2str(uniqueruns(i))];
    writetable(T,excel_filepath,'Sheet',sheet_name,'Range','A1')
end

plot_upstream_PB    = 0;
plot_upstream_SRA   = 0;
plot_downstream     = 0;
plot_repeatability  = 1;

%% Plotting Mean Velocities
unc_linewidth = 1;

if plot_repeatability
    %PB_repeatability
            labels_plot = strings(length(PB_repeat),1);
            maxheight = 0;
            figure(3);
            for j = 1:length(PB_repeat)
                i = runs_list(uniqueruns==PB_repeat(j));
                plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(PB_repeat(j)));
            end
                plot(aw.u(:,3),aw.h(:,3),cfd_colorlist(1),'Linewidth',2);
                labels_plot(4) = strcat("RANS CFD at matching Location");
            for j = 1:length(PB_repeat)
                i = runs_list(uniqueruns==PB_repeat(j));
                color_plot = [':',colorlist(j)];
                plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                maxheight = max([maxheight,max(c_heights{i})]);
            end
        
            grid on;    
            title('Pizza Box Repeatability');
            xlabel('Mean Velocity [m/s]');
            ylabel('Height above the surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([0,1000]);
            ylim([0,maxheight+1])

     %SRA_repeatability
            labels_plot = strings(length(SRA_repeat),1);
            maxheight = 0;
            figure(4);
            for j = 1:length(SRA_repeat)
                i = runs_list(uniqueruns==SRA_repeat(j));
                plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(SRA_repeat(j)));
                maxheight = max([maxheight,max(c_heights{i})]);
            end
                plot(aw.u(:,3),aw.h(:,3),cfd_colorlist(1),'Linewidth',2);
                labels_plot(4) = strcat("RANS CFD at matching Location");
            for j = 1:length(SRA_repeat)
                i = runs_list(uniqueruns==SRA_repeat(j));
                color_plot = [':',colorlist(j)];
                plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            end
        
            grid on;    
            title('Synthetic Roughness Repeatability');
            xlabel('Velocity [m/s]');
            ylabel('Height above the surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([0,1000]);
            ylim([0,maxheight+1])

end

if plot_upstream_PB
     %PB_BL seperate plots
            maxheight = 0;  
            for j= 1:length(PB_BL)
                i = runs_list(uniqueruns==PB_BL(j));
                maxheight = max([maxheight,max(c_heights{i})]);
            end

        j = 1;
        i = runs_list(uniqueruns==PB_BL(j));
        color_plot = [':',colorlist(j)];
        downstream_loc_run = downstream_loc(i);

        figure(6);
        subplot(1,3,j);
        plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(aw.u(:,j),aw.h(:,j),cfd_colorlist(1),'Linewidth',2);
                labels_plot= [strcat("FLEET-Location ",num2str(locations_list(i))),"","","RANS CFD"];
        grid on;    
%         title(['Location ',num2str(locations_list(i))]);
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,1000]);
        ylim([0,maxheight+1])
        hold off;

        j = 2;
        i = runs_list(uniqueruns==PB_BL(j));
        color_plot = [':',colorlist(j)];
        downstream_loc_run = downstream_loc(i);

        subplot(1,3,j);
        plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(aw.u(:,j),aw.h(:,j),cfd_colorlist(1),'Linewidth',2);
                labels_plot= [strcat("FLEET-Location ",num2str(locations_list(i))),"","","RANS CFD"];
        grid on;    
        title(['Location ',num2str(locations_list(i))]);
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,1000]);
        ylim([0,maxheight+1])
        hold off;

        subplot(1,3,3)
        labels_plot = strings(length(PB_repeat),1);
            for j = 1:length(PB_repeat)
                i = runs_list(uniqueruns==PB_repeat(j));
                plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("FLEET-Location ",num2str(locations_list(i)));
            end
                plot(aw.u(:,3),aw.h(:,3),cfd_colorlist(1),'Linewidth',2);
                labels_plot(4) = strcat("RANS CFD");
            for j = 1:length(PB_repeat)
                i = runs_list(uniqueruns==PB_repeat(j));
                color_plot = [':',colorlist(j)];
                plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            end
            
            downstream_loc_run = downstream_loc(i);
            grid on;    
            title(['Location ',num2str(locations_list(i))]);
            xlabel('Velocity [m/s]');
            ylabel('Height above the surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([0,1000]);
            ylim([0,maxheight+1])

    %PB_BL_sp (single plot)

        labels_plot = strings(length(PB_BL_sp),1);
        maxheight = 0;
        figure(7);
        for j = 1:length(PB_BL_sp)
            i = runs_list(uniqueruns==PB_BL_sp(j));
            plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
            hold on;
            labels_plot(j) = strcat("Location",num2str(locations_list(i)));
            maxheight = max([maxheight,max(c_heights{i})]);
        end
        for j = 1:length(PB_BL_sp)
            i = runs_list(uniqueruns==PB_BL_sp(j));
            color_plot = [':',colorlist(j)];
            plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        end
    
        grid on;    
%         title('PB Boundary Layer Evolution');
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,1000]);
        ylim([0,maxheight+1]);
end

if plot_upstream_SRA
    %SRA_BL
            maxheight = 0;  
            for j= 1:length(SRA_BL)
                i = runs_list(uniqueruns==SRA_BL(j));
                maxheight = max([maxheight,max(c_heights{i})]);
            end

        j = 1;
        i = runs_list(uniqueruns==SRA_BL(j));
        color_plot = [':',colorlist(j)];
        downstream_loc_run = downstream_loc(i);

        figure(8);
        subplot(1,3,j);
        plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        grid on;    
        title(['SRA ',num2str(downstream_loc_run),' mm']);
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,1000]);
        ylim([0,maxheight+1]);
        i_old = i;

        j = 2;
        i = runs_list(uniqueruns==SRA_BL(j));
        color_plot = [':',colorlist(j)];
        downstream_loc_run = downstream_loc(i);

        subplot(1,3,1);
        plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(aw.u(:,1),aw.h(:,1),cfd_colorlist(1),'Linewidth',2);
                labels_plot= [strcat("FLEET Run: ",num2str(uniqueruns(i_old))),"","",strcat("FLEET Run: ",num2str(uniqueruns(i))),"","","RANS CFD at matching Location"];
        grid on;    
        title(['SRA ',num2str(downstream_loc_run),' mm']);
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,1000]);
        ylim([0,maxheight+1]);
        hold off;

        j = 3;
        i = runs_list(uniqueruns==SRA_BL(j));
        color_plot = [':',colorlist(j)];
        downstream_loc_run = downstream_loc(i);

        subplot(1,3,2);
        plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(aw.u(:,2),aw.h(:,2),cfd_colorlist(1),'Linewidth',2);
                labels_plot= [strcat("FLEET Run: ",num2str(uniqueruns(i))),"","","RANS CFD at matching Location"];
        grid on;    
        title(['SRA ',num2str(downstream_loc_run),' mm']);
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,1000]);
        ylim([0,maxheight+1]);
        hold off;

        subplot(1,3,3)
        labels_plot = strings(length(SRA_repeat),1);
            for j = 1:length(SRA_repeat)
                i = runs_list(uniqueruns==SRA_repeat(j));
                plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(SRA_repeat(j)));
            end
                plot(aw.u(:,3),aw.h(:,3),cfd_colorlist(1),'Linewidth',2);
                labels_plot(4)= "RANS CFD at matching Location";
            for j = 1:length(SRA_repeat)
                i = runs_list(uniqueruns==SRA_repeat(j));
                color_plot = [':',colorlist(j)];
                plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            end
            
            downstream_loc_run = downstream_loc(i);
            grid on;    
            title(['SRA ',num2str(downstream_loc_run),' mm']);
            xlabel('Velocity [m/s]');
            ylabel('Height above the surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([0,1000]);
            ylim([0,maxheight+1]);

    %SRA_BL_sp (single plot)

            labels_plot = strings(length(SRA_BL_sp),1);
            maxheight = 0;
            figure(9);
            for j = 1:length(SRA_BL_sp)
                i = runs_list(uniqueruns==SRA_BL_sp(j));
                plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat(num2str(downstream_loc(i))," mm");
                maxheight = max([maxheight,max(c_heights{i})]);
            end
            for j = 1:length(SRA_BL_sp)
                i = runs_list(uniqueruns==SRA_BL_sp(j));
                color_plot = [':',colorlist(j)];
                plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            end
        
            grid on;    
            title('SRA Boundary Layer Evolution');
            xlabel('Velocity [m/s]');
            ylabel('Height above the surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([0,1000]);
            ylim([0,maxheight+1]);
end

if plot_downstream
    %PB_downstm seperate subplots
        
        cfd_j = [4,6,5];
        figure(10);
        maxheight = 0;
        labels_plot = strings(length(PB_downstm),1);
            for j = 1:length(PB_downstm)
                subplot(1,3,j);
                i = runs_list(uniqueruns==PB_downstm(j));
                plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(PB_downstm(j)));
                maxheight = max([maxheight,max(c_heights{i})]);
            end
            for j = 1:length(PB_downstm)
                subplot(1,3,j);
                i = runs_list(uniqueruns==PB_downstm(j));
                color_plot = [':',colorlist(j)];
                plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                    plot(aw.u(:,cfd_j(j)),aw.h(:,cfd_j(j)),cfd_colorlist(1),'Linewidth',2);
                downstream_loc_run = downstream_loc(i);
                spanwise_loc_run = spanwise_loc(i);
                grid on;    
                title(['PB ',num2str(downstream_loc_run),' mm downstream, ',num2str(spanwise_loc_run),' mm from center']);
                legend([strcat("FLEET Run: ",num2str(uniqueruns(i))),"","","RANS CFD at matching Location"])
                xlabel('Velocity [m/s]');
                ylabel('Height above the surface [mm]');
                set(gca,'FontSize', 15);
                set(gca,'fontname','times')  % Set it to times
                xlim([0,900]);
                ylim([0,maxheight+1]);
            end

    %PB_downstm single plot
        labels_plot = strings(length(PB_downstm),1);
        maxheight = 0;
        figure(11);
        for j = 1:length(PB_downstm)
            i = runs_list(uniqueruns==PB_downstm(j));
            plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
            hold on;
            downstream_loc_run = downstream_loc(i);
            spanwise_loc_run = spanwise_loc(i);
            labels_plot(j) = strcat([num2str(downstream_loc_run),' mm downstream, ',num2str(spanwise_loc_run),' mm from center']);
        end
        for j = 1:length(PB_downstm)
            i = runs_list(uniqueruns==PB_downstm(j));
            color_plot = [':',colorlist(j)];
            plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            maxheight = max([maxheight,max(c_heights{i})]);
        end
    
        grid on;    
        title('PB Downstream BL');
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,900]);
        ylim([0,maxheight+1])

end

%% Plotting RMS Velocities
unc_linewidth = 1;

if plot_repeatability

    %PB_repeatability
            labels_plot = strings(length(PB_repeat),1);
            maxheight = 0;
            figure(12);
            for j = 1:length(PB_repeat)
                i = runs_list(uniqueruns==PB_repeat(j));
                plot(c_rms{i},c_heights{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(PB_repeat(j)));
            end
            for j = 1:length(PB_repeat)
                i = runs_list(uniqueruns==PB_repeat(j));
                color_plot = [':',colorlist(j)];
                plot(c_rms{i}-unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_rms{i}+unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                maxheight = max([maxheight,max(c_heights{i})]);
            end
        
            grid on;    
            title('Pizza Box Repeatability');
            xlabel('RMS Velocity [m/s]');
            ylabel('Height above the surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([0,200]);
            ylim([0,maxheight+1])

     %SRA_repeatability
            labels_plot = strings(length(SRA_repeat),1);
            maxheight = 0;
            figure(13);
            for j = 1:length(SRA_repeat)
                i = runs_list(uniqueruns==SRA_repeat(j));
                plot(c_rms{i},c_heights{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(SRA_repeat(j)));
                maxheight = max([maxheight,max(c_heights{i})]);
            end
            for j = 1:length(SRA_repeat)
                i = runs_list(uniqueruns==SRA_repeat(j));
                color_plot = [':',colorlist(j)];
                plot(c_rms{i}-unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_rms{i}+unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            end
        
            grid on;    
            title('Synthetic Roughness Repeatability');
            xlabel('RMS Velocity [m/s]');
            ylabel('Height above the surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([0,200]);
            ylim([0,maxheight+1])

end

    %SRA vs PB at 22C
            both_repeat = [PB_repeat,SRA_repeat];
            labels_plot = ["PB Run","SRA Run","","","",""];
            proc_list = [1,4,2,3,5,6];
            maxheight = 0;
            figure(14);
            for j = proc_list
                c = ceil(j/length(PB_repeat));
                i = runs_list(uniqueruns==both_repeat(j));
                plot(c_rms{i},c_heights{i},colorlist(c),'Linewidth',2);
                hold on;
                maxheight = max([maxheight,max(c_heights{i})]);
            end
            for j = 1:length(both_repeat)
                i = runs_list(uniqueruns==both_repeat(j));
                c = ceil(j/length(PB_repeat));
                color_plot = [':',colorlist(c)];
                plot(c_rms{i}-unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_rms{i}+unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            end
        
            grid on;    
            title('Pizza Box vs Synthetic Roughness');
            xlabel('RMSVelocity [m/s]');
            ylabel('Height above the surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([0,200]);
            ylim([0,maxheight+1])

if plot_upstream_PB
    %PB_BL seperate plots
            maxheight = 0;  
            for j= 1:length(PB_BL)
                i = runs_list(uniqueruns==PB_BL(j));
                maxheight = max([maxheight,max(c_heights{i})]);
            end

        j = 1;
        i = runs_list(uniqueruns==PB_BL(j));
        color_plot = [':',colorlist(j)];
        downstream_loc_run = downstream_loc(i);

        figure(15);
        subplot(1,3,j);
        plot(c_rms{i},c_heights{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_rms{i}-unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_rms{i}+unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        labels_plot= [strcat("FLEET Run: ",num2str(uniqueruns(i)))];
        grid on;    
        title(['PB ',num2str(downstream_loc_run),' mm']);
        xlabel('RMS Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,200]);
        ylim([0,maxheight+1])
        hold off;

        j = 2;
        i = runs_list(uniqueruns==PB_BL(j));
        color_plot = [':',colorlist(j)];
        downstream_loc_run = downstream_loc(i);

        subplot(1,3,j);
        plot(c_rms{i},c_heights{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_rms{i}-unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_rms{i}+unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        labels_plot= [strcat("Location",num2str(locations_list(i)))];
        grid on;    
        title(['PB ',num2str(downstream_loc_run),' mm']);
        xlabel('RMS Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,200]);
        ylim([0,maxheight+1])
        hold off;

        subplot(1,3,3)
        labels_plot = strings(length(PB_repeat),1);
            for j = 1:length(PB_repeat)
                i = runs_list(uniqueruns==PB_repeat(j));
                plot(c_rms{i},c_heights{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(PB_repeat(j)));
            end

            for j = 1:length(PB_repeat)
                i = runs_list(uniqueruns==PB_repeat(j));
                color_plot = [':',colorlist(j)];
                plot(c_rms{i}-unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_rms{i}+unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            end
            
            downstream_loc_run = downstream_loc(i);
            grid on;    
            title(['PB ',num2str(downstream_loc_run),' mm']);
            xlabel('RMS Velocity [m/s]');
            ylabel('Height above the surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([0,200]);
            ylim([0,maxheight+1])

    %PB_BL_sp (single plot)

        labels_plot = strings(length(PB_BL_sp),1);
        maxheight = 0;
        figure(16);
        for j = 1:length(PB_BL_sp)
            i = runs_list(uniqueruns==PB_BL_sp(j));
            plot(c_rms{i},c_heights{i},colorlist(j),'Linewidth',2);
            hold on;
            labels_plot(j) = strcat("Location",num2str(locations_list(i)));
            maxheight = max([maxheight,max(c_heights{i})]);
        end
        for j = 1:length(PB_BL_sp)
            i = runs_list(uniqueruns==PB_BL_sp(j));
            color_plot = [':',colorlist(j)];
            plot(c_rms{i}-unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            plot(c_rms{i}+unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        end
    
        grid on;    
        title('PB Boundary Layer Evolution');
        xlabel('RMS Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,200]);
        ylim([0,maxheight+1]);
end

if plot_upstream_SRA

    %SRA_BL
            maxheight = 0;  
            for j= 1:length(SRA_BL)
                i = runs_list(uniqueruns==SRA_BL(j));
                maxheight = max([maxheight,max(c_heights{i})]);
            end

        j = 1;
        i = runs_list(uniqueruns==SRA_BL(j));
        color_plot = [':',colorlist(j)];
        downstream_loc_run = downstream_loc(i);

        figure(17);
        subplot(1,3,j);
        plot(c_rms{i},c_heights{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_rms{i}-unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_rms{i}+unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        grid on;    
        title(['SRA ',num2str(downstream_loc_run),' mm']);
        xlabel('RMS Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,200]);
        ylim([0,maxheight+1]);
        i_old = i;

        j = 2;
        i = runs_list(uniqueruns==SRA_BL(j));
        color_plot = [':',colorlist(j)];
        downstream_loc_run = downstream_loc(i);

        subplot(1,3,1);
        plot(c_rms{i},c_heights{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_rms{i}-unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_rms{i}+unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
               labels_plot= [strcat("FLEET Run: ",num2str(uniqueruns(i_old))),"","",strcat("FLEET Run: ",num2str(uniqueruns(i)))];
        grid on;    
        title(['SRA ',num2str(downstream_loc_run),' mm']);
        xlabel('RMS Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,200]);
        ylim([0,maxheight+1]);
        hold off;

        j = 3;
        i = runs_list(uniqueruns==SRA_BL(j));
        color_plot = [':',colorlist(j)];
        downstream_loc_run = downstream_loc(i);

        subplot(1,3,2);
        plot(c_rms{i},c_heights{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_rms{i}-unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_rms{i}+unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
               labels_plot= [strcat("FLEET Run: ",num2str(uniqueruns(i)))];
        grid on;    
        title(['SRA ',num2str(downstream_loc_run),' mm']);
        xlabel('RMS Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,200]);
        ylim([0,maxheight+1]);
        hold off;

        subplot(1,3,3)
        labels_plot = strings(length(SRA_repeat),1);
            for j = 1:length(SRA_repeat)
                i = runs_list(uniqueruns==SRA_repeat(j));
                plot(c_rms{i},c_heights{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(SRA_repeat(j)));
            end

            for j = 1:length(SRA_repeat)
                i = runs_list(uniqueruns==SRA_repeat(j));
                color_plot = [':',colorlist(j)];
                plot(c_rms{i}-unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_rms{i}+unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            end
            
            downstream_loc_run = downstream_loc(i);
            grid on;    
            title(['SRA ',num2str(downstream_loc_run),' mm']);
            xlabel('RMS Velocity [m/s]');
            ylabel('Height above the surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([0,200]);
            ylim([0,maxheight+1]);

    %SRA_BL_sp (single plot)

            labels_plot = strings(length(SRA_BL_sp),1);
            maxheight = 0;
            figure(18);
            for j = 1:length(SRA_BL_sp)
                i = runs_list(uniqueruns==SRA_BL_sp(j));
                plot(c_rms{i},c_heights{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat(num2str(downstream_loc(i))," mm");
                maxheight = max([maxheight,max(c_heights{i})]);
            end
            for j = 1:length(SRA_BL_sp)
                i = runs_list(uniqueruns==SRA_BL_sp(j));
                color_plot = [':',colorlist(j)];
                plot(c_rms{i}-unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_rms{i}+unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            end
        
            grid on;    
            title('SRA Boundary Layer Evolution');
            xlabel('RMS Velocity [m/s]');
            ylabel('Height above the surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([0,200]);
            ylim([0,maxheight+1]);

end

if plot_downstream

    %PB_downstm seperate subplots
        
        cfd_j = [4,6,5];
        figure(19);
        maxheight = 0;
        labels_plot = strings(length(PB_downstm),1);
            for j = 1:length(PB_downstm)
                subplot(1,3,j);
                i = runs_list(uniqueruns==PB_downstm(j));
                plot(c_rms{i},c_heights{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(PB_downstm(j)));
                maxheight = max([maxheight,max(c_heights{i})]);
            end
            for j = 1:length(PB_downstm)
                subplot(1,3,j);
                i = runs_list(uniqueruns==PB_downstm(j));
                color_plot = [':',colorlist(j)];
                plot(c_rms{i}-unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_rms{i}+unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                downstream_loc_run = downstream_loc(i);
                spanwise_loc_run = spanwise_loc(i);
                grid on;    
                title(['PB ',num2str(downstream_loc_run),' mm downstream, ',num2str(spanwise_loc_run),' mm from center']);
                legend([strcat("FLEET Run: ",num2str(uniqueruns(i)))])
                xlabel('RMS Velocity [m/s]');
                ylabel('Height above the surface [mm]');
                set(gca,'FontSize', 15);
                set(gca,'fontname','times')  % Set it to times
                xlim([0,200]);
                ylim([0,maxheight+1]);
            end

    %PB_downstm single plot
        labels_plot = strings(length(PB_downstm),1);
        maxheight = 0;
        figure(20);
        for j = 1:length(PB_downstm)
            i = runs_list(uniqueruns==PB_downstm(j));
            plot(c_rms{i},c_heights{i},colorlist(j),'Linewidth',2);
            hold on;
            downstream_loc_run = downstream_loc(i);
            spanwise_loc_run = spanwise_loc(i);
            labels_plot(j) = strcat([num2str(downstream_loc_run),' mm downstream, ',num2str(spanwise_loc_run),' mm from center']);
        end
        for j = 1:length(PB_downstm)
            i = runs_list(uniqueruns==PB_downstm(j));
            color_plot = [':',colorlist(j)];
            plot(c_rms{i}-unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            plot(c_rms{i}+unc_rms_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            maxheight = max([maxheight,max(c_heights{i})]);
        end
    
        grid on;    
        title('PB Downstream BL');
        xlabel('RMS Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,200]);
        ylim([0,maxheight+1])

end

%% Scitech Paper Plots
close all;

% PB_repeat   = [5484,5485,5486];
% SRA_repeat  = [5479,5483,5495];
% PB_BL       = [5491,5490,PB_repeat];
% PB_BL_sp    = [5491,5490,5486];
% SRA_BL      = [5492,5493,5494,SRA_repeat];
% SRA_BL_sp   = [5493,5494,5495];
% PB_downstm  = [5487,5488,5489]; %C, R, C-downstream
% all_runs    = [5479,5483:5495];
lnw = 2;
cfd_lnw = 2;

figure(21); %men and rms velocity upstream of the shocks
v_inf = zeros(3,1);
subplot(1,3,1); %mean, 181 mm

        i = runs_list(uniqueruns==PB_BL(1));
        plot(c_velocities{i},c_heights{i},'m','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==SRA_BL(1));
        plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(2));
        plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
        plot(aw.u(:,1),aw.h(:,1),cfd_colorlist(1),'Linewidth',cfd_lnw);
            v_inf(1) = max(aw.u(:,1));

    grid on;    
    title('Location 1');
    xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
    ylabel('Height above surface [mm]');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,925]);
    ylim([0,12]);

subplot(1,3,2); %mean, 320 mm

        i = runs_list(uniqueruns==PB_BL(2));
        plot(c_velocities{i},c_heights{i},'m','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==SRA_BL(3));
        plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
        plot(aw.u(:,2),aw.h(:,2),cfd_colorlist(1),'Linewidth',cfd_lnw);
            v_inf(2) = max(aw.u(:,2));
        labels_plot= ["PB FLEET","SRA FLEET","RANS CFD"];

    legend(labels_plot)
    grid on;    
    title('Location 2');
    xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,925]);
    ylim([0,12]);

subplot(1,3,3); %mean, 370 mm

            i = runs_list(uniqueruns==PB_BL(3));
        plot(c_velocities{i},c_heights{i},'m','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==PB_BL(4));
        plot(c_velocities{i},c_heights{i},'m','Linewidth',lnw);
            i = runs_list(uniqueruns==PB_BL(5));
        plot(c_velocities{i},c_heights{i},'m','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(4));
        plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(5));
        plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(6));
        plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
        plot(aw.u(:,3),aw.h(:,3),cfd_colorlist(1),'Linewidth',cfd_lnw);
            v_inf(3) = max(aw.u(:,3));


    grid on;    
    title('Location 3');
    xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,925]);
    ylim([0,12]);

figure(35);
subplot(1,3,1); %RMS, 181 mm

            i = runs_list(uniqueruns==PB_BL(1));
        plot(c_rms{i}./v_inf(1),c_heights{i},'m','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==SRA_BL(1));
        plot(c_rms{i}./v_inf(1),c_heights{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(2));
        plot(c_rms{i}./v_inf(1),c_heights{i},'b','Linewidth',lnw);
        plot(DNS_RMS./max(DNS_V),DNS_h,'--k','Linewidth',lnw)
    grid on;    
    xlabel('$V_{RMS} / \bar{V}_{e}$','Interpreter','Latex')
    ylabel('Height above surface [mm]');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,0.2]);
    ylim([0,6]);
    title('Location 1');


subplot(1,3,2); %RMS, 320 mm

            i = runs_list(uniqueruns==PB_BL(2));
        plot(c_rms{i}./v_inf(2),c_heights{i},'m','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==SRA_BL(3));
        plot(c_rms{i}./v_inf(2),c_heights{i},'b','Linewidth',lnw);
        plot(DNS_RMS./max(DNS_V),DNS_h,'--k','Linewidth',lnw)

    grid on;    
    legend(["PB FLEET","SRA FLEET","*DNS CFD*"])
    xlabel('$V_{RMS} / \bar{V}_{e}$','Interpreter','Latex')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,0.2]);
    ylim([0,6]);
    title('Location 2');


subplot(1,3,3); %RMS, 370 mm

            i = runs_list(uniqueruns==PB_BL(3));
        plot(c_rms{i}./v_inf(3),c_heights{i},'m','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==PB_BL(4));
        plot(c_rms{i}./v_inf(3),c_heights{i},'m','Linewidth',lnw);
            i = runs_list(uniqueruns==PB_BL(5));
        plot(c_rms{i}./v_inf(3),c_heights{i},'m','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(4));
        plot(c_rms{i}./v_inf(3),c_heights{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(5));
        plot(c_rms{i}./v_inf(3),c_heights{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(6));
        plot(c_rms{i}./v_inf(3),c_heights{i},'b','Linewidth',lnw);
            xline(0.5,':k');
        plot(DNS_RMS./max(DNS_V),DNS_h,'--k','Linewidth',lnw)

    grid on;    
    xlabel('$V_{RMS} / \bar{V}_{e}$','Interpreter','Latex')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,0.2]);
    ylim([0,6]);
    title('Location 3');


% 
%         for m = 1:1
%             lnw = 2;
%             cfd_lnw = 2;
%             
%             figure(28); %mean L1
%             v_inf = zeros(3,1);
%             
%                     i = runs_list(uniqueruns==PB_BL(1));
%                     plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
%                     hold on;
%                         i = runs_list(uniqueruns==SRA_BL(1));
%                     plot(c_velocities{i},c_heights{i},'r','Linewidth',lnw);
%                         i = runs_list(uniqueruns==SRA_BL(2));
%                     plot(c_velocities{i},c_heights{i},'r','Linewidth',lnw);
%                     plot(aw.u(:,1),aw.h(:,1),cfd_colorlist(1),'Linewidth',cfd_lnw);
%                         v_inf(1) = max(aw.u(:,1));
% 
%                 labels_plot= ["PB FLEET","SRA FLEET","","RANS CFD"];
%                 legend(labels_plot)
%                 grid on;    
%                 title('Location 1');
%                 xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
%                 ylabel('Height above surface [mm]');
%                 set(gca,'FontSize', 15);
%                 set(gca,'fontname','times')  % Set it to times
%                 xlim([0,925]);
%                 ylim([0,12]);
%             
%             figure(29); %mean L2
%             
%                     i = runs_list(uniqueruns==PB_BL(2));
%                     plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
%                     hold on;
%                         i = runs_list(uniqueruns==SRA_BL(3));
%                     plot(c_velocities{i},c_heights{i},'r','Linewidth',lnw);
%                     plot(aw.u(:,2),aw.h(:,2),cfd_colorlist(1),'Linewidth',cfd_lnw);
%                         v_inf(2) = max(aw.u(:,2));
% 
%                 labels_plot= ["PB FLEET","SRA FLEET","RANS CFD"];
%                 legend(labels_plot)
%                 grid on;    
%                 title('Location 2');
%                 xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
%                 ylabel('Height above surface [mm]');
%                 set(gca,'FontSize', 15);
%                 set(gca,'fontname','times')  % Set it to times
%                 xlim([0,925]);
%                 ylim([0,12]);
%             
%             figure(30); %mean L3
%             
%                         i = runs_list(uniqueruns==PB_BL(3));
%                     plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
%                     hold on;
%                         i = runs_list(uniqueruns==PB_BL(4));
%                     plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
%                         i = runs_list(uniqueruns==PB_BL(5));
%                     plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
%                         i = runs_list(uniqueruns==SRA_BL(4));
%                     plot(c_velocities{i},c_heights{i},'r','Linewidth',lnw);
%                         i = runs_list(uniqueruns==SRA_BL(5));
%                     plot(c_velocities{i},c_heights{i},'r','Linewidth',lnw);
%                         i = runs_list(uniqueruns==SRA_BL(6));
%                     plot(c_velocities{i},c_heights{i},'r','Linewidth',lnw);
%                     plot(aw.u(:,3),aw.h(:,3),cfd_colorlist(1),'Linewidth',cfd_lnw);
%                         v_inf(3) = max(aw.u(:,3));
%             
%                 labels_plot= ["PB FLEET","","","SRA FLEET","","","RANS CFD"];
%                 legend(labels_plot)
%                 grid on;    
%                 title('Location 3');
%                 xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
%                 ylabel('Height above surface [mm]');
%                 set(gca,'FontSize', 15);
%                 set(gca,'fontname','times')  % Set it to times
%                 xlim([0,925]);
%                 ylim([0,12]);
%             
%             figure(31); %RMS L1
%                         i = runs_list(uniqueruns==PB_BL(1));
%                     plot(c_rms{i}./v_inf(1),c_heights{i},'b','Linewidth',lnw);
%                     hold on;
%                         i = runs_list(uniqueruns==SRA_BL(1));
%                     plot(c_rms{i}./v_inf(1),c_heights{i},'r','Linewidth',lnw);
%                         i = runs_list(uniqueruns==SRA_BL(2));
%                     plot(c_rms{i}./v_inf(1),c_heights{i},'r','Linewidth',lnw);
%             
%                 labels_plot= ["PB FLEET","SRA FLEET","","RANS CFD"];
%                 legend(labels_plot)
%                 grid on;    
%                 title('Location 1');
%                 xlabel('$V_{RMS} / \bar{V}_{e}$','Interpreter','Latex')
%                 ylabel('Height above surface [mm]');
%                 set(gca,'FontSize', 15);
%                 set(gca,'fontname','times')  % Set it to times
%                 xlim([0,0.2]);
%                 ylim([0,6]);
%             
%             figure(32); %RMS L2
%             
%                         i = runs_list(uniqueruns==PB_BL(2));
%                     plot(c_rms{i}./v_inf(2),c_heights{i},'b','Linewidth',lnw);
%                     hold on;
%                         i = runs_list(uniqueruns==SRA_BL(3));
%                     plot(c_rms{i}./v_inf(2),c_heights{i},'r','Linewidth',lnw);
%             
%                 labels_plot= ["PB FLEET","SRA FLEET"];
%                 legend(labels_plot)
%                 grid on;    
%                 title('Location 2');
%                 xlabel('$V_{RMS} / \bar{V}_{e}$','Interpreter','Latex')
%                 ylabel('Height above surface [mm]');
%                 set(gca,'FontSize', 15);
%                 set(gca,'fontname','times')  % Set it to times
%                 xlim([0,0.2]);
%                 ylim([0,6]);
%             
%             figure(33); %RMS L3
%             
%                         i = runs_list(uniqueruns==PB_BL(3));
%                     plot(c_rms{i}./v_inf(3),c_heights{i},'b','Linewidth',lnw);
%                     hold on;
%                         i = runs_list(uniqueruns==PB_BL(4));
%                     plot(c_rms{i}./v_inf(3),c_heights{i},'b','Linewidth',lnw);
%                         i = runs_list(uniqueruns==PB_BL(5));
%                     plot(c_rms{i}./v_inf(3),c_heights{i},'b','Linewidth',lnw);
%                         i = runs_list(uniqueruns==SRA_BL(4));
%                     plot(c_rms{i}./v_inf(3),c_heights{i},'r','Linewidth',lnw);
%                         i = runs_list(uniqueruns==SRA_BL(5));
%                     plot(c_rms{i}./v_inf(3),c_heights{i},'r','Linewidth',lnw);
%                         i = runs_list(uniqueruns==SRA_BL(6));
%                     plot(c_rms{i}./v_inf(3),c_heights{i},'r','Linewidth',lnw);
%                         xline(0.5,':k');
%             
%                 labels_plot= ["PB FLEET","","","SRA FLEET","",""];
%                 legend(labels_plot)
%                 title('Location 3');
%                 grid on;    
%                 xlabel('$V_{RMS} / \bar{V}_{e}$','Interpreter','Latex')
%                 ylabel('Height above surface [mm]');
%                 set(gca,'FontSize', 15);
%                 set(gca,'fontname','times')  % Set it to times
%                 xlim([0,0.2]);
%                 ylim([0,6]);
%         
%         
%         end
% 



    figure(22); %uncertainty upstream of the shocks
    lnw = 2;
        subplot(1,2,1);
            i = runs_list(uniqueruns==PB_BL(3));
        plot(unc_velo_COMBINED{i},c_heights{i},'b','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==PB_BL(4));
        plot(unc_velo_COMBINED{i},c_heights{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==PB_BL(5));
        plot(unc_velo_COMBINED{i},c_heights{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(4));
        plot(unc_velo_COMBINED{i},c_heights{i},'r','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(5));
        plot(unc_velo_COMBINED{i},c_heights{i},'r','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(6));
        plot(unc_velo_COMBINED{i},c_heights{i},'r','Linewidth',lnw);
        yline(0.4,':k','Linewidth',2);
        labels_plot= ["PB FLEET","","","SRA FLEET","","","Refl. Light Height"];

    grid on;    
    xlabel('Uncertainty in $\bar{V}$ [m/s]','Interpreter','Latex')
    ylabel('Height above surface [mm]');
    legend(labels_plot)
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,150]);
    ylim([0,8]);

      subplot(1,2,2);
            i = runs_list(uniqueruns==PB_BL(3));
        plot(unc_rms_COMBINED{i}./v_inf(3),c_heights{i},'b','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==PB_BL(4));
        plot(unc_rms_COMBINED{i}./v_inf(3),c_heights{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==PB_BL(5));
        plot(unc_rms_COMBINED{i}./v_inf(3),c_heights{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(4));
        plot(unc_rms_COMBINED{i}./v_inf(3),c_heights{i},'r','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(5));
        plot(unc_rms_COMBINED{i}./v_inf(3),c_heights{i},'r','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(6));
        plot(unc_rms_COMBINED{i}./v_inf(3),c_heights{i},'r','Linewidth',lnw);
        yline(0.4,':k','Linewidth',2);

    grid on;    
    xlabel('Uncertainty in $V_{RMS} / \bar{V}_{e}$','Interpreter','Latex')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,0.1]);
    ylim([0,8]);

    figure(23); %downstream velocity and RMS velocity
        cfd_j = [4,6,5];
        maxheight = 0;
        labels_plot = strings(length(PB_downstm),2);
            for j = 1:length(PB_downstm)
                subplot(1,3,j);
                i = runs_list(uniqueruns==PB_downstm(j));
                plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
                hold on;
                maxheight = max([maxheight,max(c_heights{i})]);
            end
            for j = 1:length(PB_downstm)
                subplot(1,3,j);
                i = runs_list(uniqueruns==PB_downstm(j));
                color_plot = ['--',colorlist(j)];
                plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                    plot(aw.u(:,cfd_j(j)),aw.h(:,cfd_j(j)),cfd_colorlist(1),'Linewidth',2);
                downstream_loc_run = downstream_loc(i);
                spanwise_loc_run = spanwise_loc(i);
                grid on;    
                legend('PB FLEET','','','RANS CFD')
                xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
                ylabel('Height above surface [mm]');
                set(gca,'FontSize', 15);
                set(gca,'fontname','times')  % Set it to times
                xlim([0,900]);
                ylim([0,10]);
                title(strcat("Location ",num2str(locations_list(j+5))))
            end

    mean_upstream_vinf = 862.457;
    figure(24);
      labels_plot = strings(length(PB_downstm),1);
        maxheight = 0;
        for j = 1:length(PB_downstm)
            i = runs_list(uniqueruns==PB_downstm(j));
            plot(c_rms{i}./mean_upstream_vinf,c_heights{i},colorlist(j),'Linewidth',2);
            hold on;
            downstream_loc_run = downstream_loc(i);
            spanwise_loc_run = spanwise_loc(i);
            labels_plot(j) = strcat("Location ",num2str(j+3));
        end
        for j = 1:length(PB_downstm)
            i = runs_list(uniqueruns==PB_downstm(j));
            color_plot = ['--',colorlist(j)];
            plot((c_rms{i}-unc_velo_COMBINED{i})./mean_upstream_vinf,c_heights{i},color_plot,'Linewidth',unc_linewidth);
            plot((c_rms{i}+unc_velo_COMBINED{i})./mean_upstream_vinf,c_heights{i},color_plot,'Linewidth',unc_linewidth);
            maxheight = max([maxheight,max(c_heights{i})]);
        end
    
        grid on;    
        xlabel('$V_{RMS} / \bar{V}_{e}^*$','Interpreter','Latex')
        ylabel('Height above surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,0.2]);
        ylim([0,7])


    %SRA vs PB at Location 3
            both_repeat = [PB_repeat(:);SRA_repeat(:)]; %3,1
            maxheight = 0;
            figure(5);
            for j = 1:length(both_repeat)
                if  ismember(both_repeat(j),PB_repeat) %pizza box
                    c = 2;
                else
                    c = 1;
                end
                i = runs_list(uniqueruns==both_repeat(j));
                plot(c_velocities{i},c_heights{i},colorlist(c),'Linewidth',2);
                hold on;
                maxheight = max([maxheight,max(c_heights{i})]);
            end
                plot(aw.u(:,3),aw.h(:,3),cfd_colorlist(1),'Linewidth',2);
            for j = 1:length(both_repeat)
                i = runs_list(uniqueruns==both_repeat(j));
                if  ismember(both_repeat(j),PB_repeat) %pizza box
                    c = 2;
                else
                    c = 1;
                end
                color_plot = ['--',colorlist(c)];
                plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            end
            plot(aw.u(:,3),aw.h(:,3),cfd_colorlist(1),'Linewidth',2);

            clear labels_plot
            labels_plot = ["PB","","","SRA","","",strcat("RANS CFD"),"PB 95% CI","","","","","","SRA 95% CI",""];
        
            grid on;    
%             title('Pizza Box vs Synthetic Roughness');
            xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
            ylabel('Height above surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([150,900]);
            ylim([0,10])
            hold on;

            

%SRA vs PB at Location 3
        both_repeat = [PB_repeat(:)]; %3,1
        maxheight = 0;
        figure(27);
        subplot(1,2,1);
        cols_repeat = [1, 0, 0; 0, 0, 1; 1, 0, 1];
        for j = 1:length(both_repeat)
            if  ismember(both_repeat(j),PB_repeat) %pizza box
                c = 2;
            else
                c = 1;
            end
            i = runs_list(uniqueruns==both_repeat(j));
            plot(c_velocities{i},c_heights{i},'Color', cols_repeat(j,:), 'LineStyle','-','Linewidth',3);
            hold on;
            maxheight = max([maxheight,max(c_heights{i})]);
        end
        for j = 1:length(both_repeat)
            i = runs_list(uniqueruns==both_repeat(j));
            if  ismember(both_repeat(j),PB_repeat) %pizza box
                c = 2;
            else
                c = 1;
            end
            plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},'Color', cols_repeat(j,:), 'LineStyle','--','Linewidth',2);
            plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},'Color', cols_repeat(j,:), 'LineStyle','--','Linewidth',2);
        end
        plot(aw.u(:,3),aw.h(:,3),cfd_colorlist(1),'Linewidth',2);

        clear labels_plot
        labels_plot = ["PB FLEET","","","95% CI","","","","","","RANS CFD"];
    
        grid on;    
%             title('Pizza Box vs Synthetic Roughness');
        xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
        ylabel('Height above surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 20);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,900]);
        ylim([0,10])
        hold on;

        subplot(1,2,2);
        both_repeat = [SRA_repeat(:)]; %3,1
        maxheight = 0;
        for j = 1:length(both_repeat)
            if  ismember(both_repeat(j),PB_repeat) %pizza box
                c = 2;
            else
                c = 1;
            end
            i = runs_list(uniqueruns==both_repeat(j));
            plot(c_velocities{i},c_heights{i},'Color', cols_repeat(j,:), 'LineStyle','-','Linewidth',3);
            hold on;
            maxheight = max([maxheight,max(c_heights{i})]);
        end
        for j = 1:length(both_repeat)
            i = runs_list(uniqueruns==both_repeat(j));
            if  ismember(both_repeat(j),PB_repeat) %pizza box
                c = 2;
            else
                c = 1;
            end
            color_plot = ['--',colorlist(c)];
            plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},'Color', cols_repeat(j,:), 'LineStyle','--','Linewidth',2);
            plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},'Color', cols_repeat(j,:), 'LineStyle','--','Linewidth',2);
        end
        plot(aw.u(:,3),aw.h(:,3),cfd_colorlist(1),'Linewidth',2);

        clear labels_plot
        labels_plot = ["SRA FLEET","","","95% CI","","","","","","RANS CFD"];
    
        grid on;    
%             title('Pizza Box vs Synthetic Roughness');
        xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
        ylabel('Height above surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 20);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,900]);
        ylim([0,10])
        hold on;

%SRA vs PB at Location 3 (RMS)
        both_repeat = [PB_repeat(:)]; %3,1
        maxheight = 0;
        figure(32);
        im_height = 10;
        subplot(1,2,1);
        for j = 1:length(both_repeat)
            if  ismember(both_repeat(j),PB_repeat) %pizza box
                c = 2;
            else
                c = 1;
            end
            i = runs_list(uniqueruns==both_repeat(j));
            plot(c_rms{i}./v_inf(3),c_heights{i},'Color', cols_repeat(j,:), 'LineStyle','-','Linewidth',3);
            hold on;
            maxheight = max([maxheight,max(c_heights{i})]);
        end
        for j = 1:length(both_repeat)
            i = runs_list(uniqueruns==both_repeat(j));
            if  ismember(both_repeat(j),PB_repeat) %pizza box
                c = 2;
            else
                c = 1;
            end
            color_plot = ['--',colorlist(c)];
            plot(c_rms{i}./v_inf(3)+unc_rms_COMBINED{i}./v_inf(3),c_heights{i},'Color', cols_repeat(j,:), 'LineStyle','--','Linewidth',2);
            plot(c_rms{i}./v_inf(3)-unc_rms_COMBINED{i}./v_inf(3),c_heights{i},'Color', cols_repeat(j,:), 'LineStyle','--','Linewidth',2);
        end
        plot(DNS_RMS./max(DNS_V),DNS_h,'-k','Linewidth',2)

        clear labels_plot
        labels_plot = ["PB FLEET","","","95% CI","","","","","","*DNS CFD*"];
    
        grid on;    
%             title('Pizza Box vs Synthetic Roughness');
        xlabel('$V_{RMS}/V_e$','Interpreter','Latex')
        ylabel('Height above surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 20);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,.25]);
        ylim([0,im_height])
        hold on;

         both_repeat = [SRA_repeat(:)]; %3,1
        maxheight = 0;
        subplot(1,2,2);
        for j = 1:length(both_repeat)
            if  ismember(both_repeat(j),PB_repeat) %pizza box
                c = 2;
            else
                c = 1;
            end
            i = runs_list(uniqueruns==both_repeat(j));
            plot(c_rms{i}./v_inf(3),c_heights{i},'Color', cols_repeat(j,:), 'LineStyle','-','Linewidth',3);
            hold on;
            maxheight = max([maxheight,max(c_heights{i})]);
        end
        for j = 1:length(both_repeat)
            i = runs_list(uniqueruns==both_repeat(j));
            if  ismember(both_repeat(j),PB_repeat) %pizza box
                c = 2;
            else
                c = 1;
            end
            color_plot = ['--',colorlist(c)];
            plot(c_rms{i}./v_inf(3)+unc_rms_COMBINED{i}./v_inf(3),c_heights{i},'Color', cols_repeat(j,:), 'LineStyle','--','Linewidth',2);
            plot(c_rms{i}./v_inf(3)-unc_rms_COMBINED{i}./v_inf(3),c_heights{i},'Color', cols_repeat(j,:), 'LineStyle','--','Linewidth',2);
        end
        plot(DNS_RMS./max(DNS_V),DNS_h,'-k','Linewidth',2)

        clear labels_plot
        labels_plot = ["SRA FLEET","","","95% CI","","","","","","*DNS CFD*"];
    
        grid on;    
%             title('Pizza Box vs Synthetic Roughness');
        xlabel('$V_{RMS}/V_e$','Interpreter','Latex')
        ylabel('Height above surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 20);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,0.25]);
        ylim([0,im_height])
        hold on;
% 
%     both_repeat = [PB_repeat(:);SRA_repeat(:)]; 
%     subplot(1,3,3);
%     for j = 1:length(both_repeat)
%             if  ismember(both_repeat(j),PB_repeat) %pizza box
%                 colp = 'm';
%             else
%                 colp = 'b';
%             end
%             i = runs_list(uniqueruns==both_repeat(j));
% 
%             mean_SNR_binary = (c_heights{i}<3)&(c_heights{i}>.5);
%             dat_snr = c_SNRs{i};
%             mean_im_snr = mean(dat_snr(mean_SNR_binary));
%             if i==1
%             colp = strcat('--',colp);
%             plot(c_SNRs{i}./mean_im_snr,c_heights{i},colp,'Linewidth',1.5);
%             else
%             plot(c_SNRs{i}./mean_im_snr,c_heights{i},colp,'Linewidth',1.5);
%             end
%             hold on;
%             maxheight = max([maxheight,max(c_heights{i})]);
%     end
% 
%         grid on;    
%         xlabel('$SNR/SNR_{max}$','Interpreter','Latex')
%         ylabel('Height above surface [mm]');
%         set(gca,'FontSize', 15);
%         set(gca,'fontname','times')  % Set it to times
%         xlim([0,2]);
%         ylim([0,im_height])
%         legend("PB, 500 mm Lens","","","SRA, 750 mm Lens","SRA, 500 mm Lens")
%         hold on;


% 
% %% Calculating t values for repeatability
% points = 200;
% h1 = zeros(points,3);
% v1 = zeros(points,3);
% s = zeros(points,3);
%     for j = 1:length(PB_repeat)
%         i = runs_list(uniqueruns==PB_repeat(j));
%         v = c_velocities{i};
%         h = c_heights{i};
%         h1(:,j) = linspace(0,6.62839742571926,points);
%         v1(:,j) = interp1(h,v,h1(:,j));
%     end
% 
%     t_crit = tinv(0.90,3);
%     h_mean = mean(h1,2);
%     v_mean = mean(v1,2);
% 
% for k = 1:length(PB_repeat)
%     for j = 1:length(h_mean)
%         s(j,k) = (v_mean(j)-v1(j,k))/(t_crit/sqrt(3));
%     end
% end
% 
% points = 200;
% h1 = zeros(points,3);
% v1 = zeros(points,3);
% s2 = zeros(points,3);
%     for j = 1:length(SRA_repeat)
%         i = runs_list(uniqueruns==SRA_repeat(j));
%         v = c_velocities{i};
%         h = c_heights{i};
%         h1(:,j) = linspace(0,6.62839742571926,points);
%         v1(:,j) = interp1(h,v,h1(:,j));
%     end
% 
%     t_crit = tinv(0.90,3);
%     h_mean = mean(h1,2);
%     v_mean = mean(v1,2);
% 
% for k = 1:length(SRA_repeat)
%     for j = 1:length(h_mean)
%         s2(j,k) = (v_mean(j)-v1(j,k))/(t_crit/sqrt(3));
%     end
% end
% 

%     %Plot RMS velocity at L3 and downstream of the shocks on the same
%     %plot
%         mean_upstream_vinf = 862.457;
%         figure(26);
%                         i = runs_list(uniqueruns==PB_BL(3));
%                         col = 'm';
%                     plot(c_rms{i}./v_inf(3),c_heights{i},col,'Linewidth',lnw);
%                     hold on;
%                         i = runs_list(uniqueruns==PB_BL(4));
%                     plot(c_rms{i}./v_inf(3),c_heights{i},col,'Linewidth',lnw);
%                         i = runs_list(uniqueruns==PB_BL(5));
%                     plot(c_rms{i}./v_inf(3),c_heights{i},col,'Linewidth',lnw);
%                         i = runs_list(uniqueruns==SRA_BL(4));
%                     plot(c_rms{i}./v_inf(3),c_heights{i},col,'Linewidth',lnw);
%                         i = runs_list(uniqueruns==SRA_BL(5));
%                     plot(c_rms{i}./v_inf(3),c_heights{i},col,'Linewidth',lnw);
%                         i = runs_list(uniqueruns==SRA_BL(6));
%                     plot(c_rms{i}./v_inf(3),c_heights{i},col,'Linewidth',lnw);
%                         xline(0.5,':k');
%                 grid on;    
%                 xlabel('$V_{RMS} / \bar{V}_{e}$','Interpreter','Latex')
%                 set(gca,'FontSize', 15);
%                 set(gca,'fontname','times')  % Set it to times
%                 xlim([0,0.2]);
%                 ylim([0,6]);
%     
%           labels_plot = strings(length(PB_downstm),1);
%             maxheight = 0;
%             for j = 1:length(PB_downstm)
%                 i = runs_list(uniqueruns==PB_downstm(j));
%                 plot(c_rms{i}./mean_upstream_vinf,c_heights{i},colorlist(j),'Linewidth',2);
%                 hold on;
%                 downstream_loc_run = downstream_loc(i);
%                 spanwise_loc_run = spanwise_loc(i);
%                 labels_plot(j) = strcat("Location ",num2str(j+3));
%             end
%             for j = 1:length(PB_downstm)
%                 i = runs_list(uniqueruns==PB_downstm(j));
%                 color_plot = ['--',colorlist(j)];
%                 plot((c_rms{i}-unc_velo_COMBINED{i})./mean_upstream_vinf,c_heights{i},color_plot,'Linewidth',unc_linewidth);
%                 plot((c_rms{i}+unc_velo_COMBINED{i})./mean_upstream_vinf,c_heights{i},color_plot,'Linewidth',unc_linewidth);
%                 maxheight = max([maxheight,max(c_heights{i})]);
%             end
%         
% 
%             labels_plot = ["Location 3";"";"";"";"";"";"";labels_plot];
%             grid on;    
%             xlabel('$V_{RMS} / \bar{V}_{e}^*$','Interpreter','Latex')
%             ylabel('Height above surface [mm]');
%             legend(labels_plot(:));
%             set(gca,'FontSize', 15);
%             set(gca,'fontname','times')  % Set it to times
%             xlim([0,0.2]);
%             ylim([0,6])

%% Inner variable plotting on the RANS CFD at location 3
Mach = 5.7;
T_inf = 58; %K
P_inf = 440; %Pa
U_inf = max(aw.u(:,3));
claus_tw = 13;
[u_plus_c,u_vd_plus_c,y_plus_c] = InnerVariableCalculator(aw.u(:,3),aw.h(:,3),Mach,T_inf,P_inf,U_inf,claus_tw);

%PB
[~,idx]=ismember(uniqueruns,PB_repeat);
PB_repeats = runs_list(idx>0);
[~,u_p_vd_PB1,y_p_vd_PB1,u_tau_c] = InnerVariableCalculator(c_velocities{PB_repeats(1)},c_heights{PB_repeats(1)},...
    Mach,T_inf,P_inf,U_inf,claus_tw);
[~,u_p_vd_PB2,y_p_vd_PB2,~] = InnerVariableCalculator(c_velocities{PB_repeats(2)},c_heights{PB_repeats(2)},...
    Mach,T_inf,P_inf,U_inf,claus_tw);
[~,u_p_vd_PB3,y_p_vd_PB3,~] = InnerVariableCalculator(c_velocities{PB_repeats(3)},c_heights{PB_repeats(3)},...
    Mach,T_inf,P_inf,U_inf,claus_tw);

    %bounds
        %upper
    [~,u_p_vd_PB_ub,y_p_vd_PB_ub] = InnerVariableCalculator(c_velocities{PB_repeats(3)}+unc_velo_COMBINED{PB_repeats(3)},...
        c_heights{PB_repeats(3)},Mach,T_inf,P_inf,U_inf,claus_tw);
        %lower
    [~,u_p_vd_PB_lb,y_p_vd_PB_lb] = InnerVariableCalculator(c_velocities{PB_repeats(1)}-unc_velo_COMBINED{PB_repeats(1)},...
        c_heights{PB_repeats(1)},Mach,T_inf,P_inf,U_inf,claus_tw);

%SRA
[~,idx]=ismember(uniqueruns,SRA_repeat);
SRA_repeats = runs_list(idx>0);
[~,u_p_vd_SRA1,y_p_vd_SRA1,~] = InnerVariableCalculator(c_velocities{SRA_repeats(1)},c_heights{SRA_repeats(1)},...
    Mach,T_inf,P_inf,U_inf,claus_tw);
[~,u_p_vd_SRA2,y_p_vd_SRA2,~] = InnerVariableCalculator(c_velocities{SRA_repeats(2)},c_heights{SRA_repeats(2)},...
    Mach,T_inf,P_inf,U_inf,claus_tw);
[~,u_p_vd_SRA3,y_p_vd_SRA3,~] = InnerVariableCalculator(c_velocities{SRA_repeats(3)},c_heights{SRA_repeats(3)},...
    Mach,T_inf,P_inf,U_inf,claus_tw);

    %bounds
        %upper
    [~,u_p_vd_SRA_ub,y_p_vd_SRA_ub] = InnerVariableCalculator(c_velocities{SRA_repeats(3)}+unc_velo_COMBINED{SRA_repeats(3)},...
        c_heights{SRA_repeats(3)},Mach,T_inf,P_inf,U_inf,claus_tw);
        %lower
    [~,u_p_vd_SRA_lb,y_p_vd_SRA_lb] = InnerVariableCalculator(c_velocities{SRA_repeats(2)}-unc_velo_COMBINED{SRA_repeats(2)},...
        c_heights{SRA_repeats(2)},Mach,T_inf,P_inf,U_inf,claus_tw);

        y_ideal = logspace(1,2.5);
        xi = 0.41;
        C = 5;
        u_ideal = (1/xi).*log(y_ideal)+C;
            %viscous sublayer
        y_visc = logspace(-1,1);
        u_visc = y_visc;
    
    %plotting innver variable cfd
        figure;
        semilogx(y_ideal,u_ideal,'k','Linewidth',2);
        hold on;
        semilogx(y_visc,u_visc,'--k','Linewidth',2);
        semilogx(y_plus_c,u_vd_plus_c,'r','Linewidth',2);
%         semilogx(DNS_yplus,DNS_uplus,'g','Linewidth',2);
        semilogx(y_p_vd_PB1,u_p_vd_PB1,'m','Linewidth',2);
        semilogx(y_p_vd_PB2,u_p_vd_PB2,'m','Linewidth',2);
        semilogx(y_p_vd_PB3,u_p_vd_PB3,'m','Linewidth',2);
        %bounds
        semilogx(y_p_vd_PB_ub,u_p_vd_PB_ub,'--m','Linewidth',1);
        semilogx(y_p_vd_PB_lb,u_p_vd_PB_lb,'--m','Linewidth',1); 

        legend('Logarithmic Law of the Wall (0.41,5)','Viscous Sublayer','RANS CFD','PB')
        set(gca,'FontSize', 20);
        set(gca,'fontname','times')  % Set it to times
        xlim([0.4,200])
        ylim([0,20])
        grid on;
        xlabel('$y^+$','Interpreter','Latex')
        ylabel('$u_{VD}^+$','Interpreter','Latex')

      %plotting innver variable cfd
        figure;
        semilogx(y_ideal,u_ideal,'k','Linewidth',2);
        hold on;
        semilogx(y_visc,u_visc,'--k','Linewidth',2);
        semilogx(y_plus_c,u_vd_plus_c,'r','Linewidth',2);
%         semilogx(DNS_yplus,DNS_uplus,'g','Linewidth',2);
        semilogx(y_p_vd_SRA1,u_p_vd_SRA1,'b','Linewidth',2);
        semilogx(y_p_vd_SRA2,u_p_vd_SRA2,'b','Linewidth',2);
        semilogx(y_p_vd_SRA3,u_p_vd_SRA3,'b','Linewidth',2);
        %bounds
        semilogx(y_p_vd_SRA_ub,u_p_vd_SRA_ub,'--b','Linewidth',1);
        semilogx(y_p_vd_SRA_lb,u_p_vd_SRA_lb,'--b','Linewidth',1);

        legend('Logarithmic Law of the Wall (0.41,5)','Viscous Sublayer','RANS CFD','SRA')
        set(gca,'FontSize', 20);
        set(gca,'fontname','times')  % Set it to times
        xlim([0.4,200])
        ylim([0,20])
        grid on;
        xlabel('$y^+$','Interpreter','Latex')
        ylabel('$u_{VD}^+$','Interpreter','Latex')

    %% Outer Variables

    bl_height = 9;
%CFD
[norm_velo_defect_RANS,y_norm_RANS,~,~] = ...
    OuterVariables(u_tau_c,aw.u(:,3),0,max(aw.u(:,3)),...
    aw.h(:,3),bl_height);

%PB
[norm_velo_defect_PB1,y_norm_PB1,norm_velo_defect_max_PB1,norm_velo_defect_min_PB1] = ...
    OuterVariables(u_tau_c,c_velocities{PB_repeats(1)},unc_velo_COMBINED{PB_repeats(1)},max(aw.u(:,3)),...
    c_heights{PB_repeats(1)},bl_height);
[norm_velo_defect_PB2,y_norm_PB2,norm_velo_defect_max_PB2,norm_velo_defect_min_PB2] = ...
    OuterVariables(u_tau_c,c_velocities{PB_repeats(2)},unc_velo_COMBINED{PB_repeats(2)},max(aw.u(:,3)),...
    c_heights{PB_repeats(2)},bl_height);
[norm_velo_defect_PB3,y_norm_PB3,norm_velo_defect_max_PB3,norm_velo_defect_min_PB3] = ...
    OuterVariables(u_tau_c,c_velocities{PB_repeats(3)},unc_velo_COMBINED{PB_repeats(3)},max(aw.u(:,3)),...
    c_heights{PB_repeats(3)},bl_height);

%SRA
[norm_velo_defect_SRA1,y_norm_SRA1,norm_velo_defect_max_SRA1,norm_velo_defect_min_SRA1] = ...
    OuterVariables(u_tau_c,c_velocities{SRA_repeats(1)},unc_velo_COMBINED{SRA_repeats(1)},max(aw.u(:,3)),...
    c_heights{SRA_repeats(1)},bl_height);
[norm_velo_defect_SRA2,y_norm_SRA2,norm_velo_defect_max_SRA2,norm_velo_defect_min_SRA2] = ...
    OuterVariables(u_tau_c,c_velocities{SRA_repeats(2)},unc_velo_COMBINED{SRA_repeats(2)},max(aw.u(:,3)),...
    c_heights{SRA_repeats(2)},bl_height);
[norm_velo_defect_SRA3,y_norm_SRA3,norm_velo_defect_max_SRA3,norm_velo_defect_min_SRA3] = ...
    OuterVariables(u_tau_c,c_velocities{SRA_repeats(3)},unc_velo_COMBINED{SRA_repeats(3)},max(aw.u(:,3)),...
    c_heights{SRA_repeats(3)},bl_height);
   
    %ideal outer variables (for no pressure gradient flat plate)
    k = 0.41;
    A = 0;
    y_delta_ideal = linspace(0,1,100);
    % velo_defect_ideal = 9.6.*((1-y_delta_ideal).^2);
    velo_defect_ideal = -(1/k).*log(y_delta_ideal)+A;
    y_delta_ideal = [y_delta_ideal,linspace(1,2,100)];
    velo_defect_ideal = [velo_defect_ideal,zeros(1,100)];
    
    %Plotting
    figure;
    plot(y_norm_RANS,norm_velo_defect_RANS,'k','Linewidth',2);
    hold on;
    %PB
    plot(y_norm_PB1,norm_velo_defect_PB1,'m','Linewidth',2);
    plot(y_norm_PB2,norm_velo_defect_PB2,'m','Linewidth',2);
    plot(y_norm_PB3,norm_velo_defect_PB3,'m','Linewidth',2);
    plot(y_norm_SRA1,norm_velo_defect_SRA1,'b','Linewidth',2);
    plot(y_norm_SRA2,norm_velo_defect_SRA2,'b','Linewidth',2);
    plot(y_norm_SRA3,norm_velo_defect_SRA3,'b','Linewidth',2);
%     
%     %bounds
%     plot(y_norm_PB3,norm_velo_defect_max_PB3,'--m','Linewidth',1);
%     plot(y_norm_PB1,norm_velo_defect_min_PB1,'--m','Linewidth',1);
%     plot(y_norm_SRA3,norm_velo_defect_max_SRA3,'--b','Linewidth',1);
%     plot(y_norm_SRA2,norm_velo_defect_min_SRA2,'--b','Linewidth',1);

    xlabel('$y/\delta$','Interpreter','Latex')
    ylabel('$\frac{(U_e-\bar{u})}{v^*}$','Interpreter','Latex')
    legend(["RANS CFD","Pizza Box","","","Synthetic Roughness"])
    set(gca,'FontSize', 18);
    set(gca,'fontname','times')  % Set it to times
    % title('Outer Variables for Turbulent Flow')
    grid on;
    ylim([-2,15])
    xlim([0,1.25])
    set(gca, 'YDir','reverse')

%% y/delta plot with snr
[ d, max_height_ind ] = min( abs( y_norm_RANS-1.2 ) );
figure;
plot(aw.u(:,3),y_norm_RANS,'k','Linewidth',2);
hold on;
%PB
plot(c_velocities{PB_repeats(1)},y_norm_PB1,'m','Linewidth',2);
plot(c_velocities{PB_repeats(2)},y_norm_PB2,'m','Linewidth',2);
plot(c_velocities{PB_repeats(3)},y_norm_PB3,'m','Linewidth',2);
plot(c_velocities{SRA_repeats(1)},y_norm_SRA1,'b','Linewidth',2);
plot(c_velocities{SRA_repeats(2)},y_norm_SRA2,'b','Linewidth',2);
plot(c_velocities{SRA_repeats(3)},y_norm_SRA3,'b','Linewidth',2);
ylabel('$y/\delta$','Interpreter','Latex')
xlabel('$\bar{V}$','Interpreter','Latex')
ylim([0,y_norm_RANS(max_height_ind)]);

yyaxis right
ylabel('$y^+$','Interpreter','Latex')
plot(aw.u(:,3),y_plus_c,'k','Linewidth',2);
ylim([0,y_plus_c(max_height_ind)]);
ax = gca;
ax.YColor = 'k';
set(gca,'FontSize', 18);
set(gca,'fontname','times')  % Set it to times
grid on;
legend(["RANS CFD","PB FLEET","","","SRA FLEET"])


%SNR

y_delta_norm = cell(length(uniqueruns),1);

i = runs_list(uniqueruns==both_repeat(1));
y_delta_norm{i} = y_norm_PB1;

i = runs_list(uniqueruns==both_repeat(2));
y_delta_norm{i} = y_norm_PB2;

i = runs_list(uniqueruns==both_repeat(3));
y_delta_norm{i} = y_norm_PB3;

i = runs_list(uniqueruns==both_repeat(4));
y_delta_norm{i} = y_norm_SRA1;

i = runs_list(uniqueruns==both_repeat(5));
y_delta_norm{i} = y_norm_SRA2;

i = runs_list(uniqueruns==both_repeat(6));
y_delta_norm{i} = y_norm_SRA3;

figure;
    for j = 1:length(both_repeat)
            if  ismember(both_repeat(j),PB_repeat) %pizza box
                colp = 'm';
            else
                colp = 'b';
            end
            i = runs_list(uniqueruns==both_repeat(j));

            if i==1
            colp = strcat('--',colp);
            plot(c_SNRs{i},y_delta_norm{i},colp,'Linewidth',2);
            else
            plot(c_SNRs{i},y_delta_norm{i},colp,'Linewidth',2);
            end
            hold on;
    end

        grid on;    
        xlabel('$SNR$','Interpreter','Latex')
        ylabel('$y/ \delta$','Interpreter','Latex')
        set(gca,'FontSize', 18);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,100]);
        ylim([0,1.2])
        legend("PB, 500 mm Lens","","","SRA, 750 mm Lens","SRA, 500 mm Lens")
        hold on;

%uncertainty
  figure; %repeatability uncertainty
    lnw = 2;
            i = runs_list(uniqueruns==PB_repeat(1));
        plot(unc_velo_COMBINED{i},y_delta_norm{i},'m','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==PB_repeat(2));
        plot(unc_velo_COMBINED{i},y_delta_norm{i},'m','Linewidth',lnw);
            i = runs_list(uniqueruns==PB_repeat(3));
        plot(unc_velo_COMBINED{i},y_delta_norm{i},'m','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_repeat(1));
        plot(unc_velo_COMBINED{i},y_delta_norm{i},'--b','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_repeat(2));
        plot(unc_velo_COMBINED{i},y_delta_norm{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_repeat(3));
        plot(unc_velo_COMBINED{i},y_delta_norm{i},'b','Linewidth',lnw);
        yline(0.5/bl_height,'k','Linewidth',2);
        labels_plot= ["PB FLEET","","","SRA 750 mm Lens","SRA 500 mm Lens","","0.5 mm"];

    grid on;    
    xlabel('Uncertainty in $\bar{V}$ [m/s]','Interpreter','Latex')
    ylabel('$y/ \delta$','Interpreter','Latex')
    legend(labels_plot)
    set(gca,'FontSize', 18);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,150]);
    ylim([0,1.2]);



%% Comparing flat plate and inlet wall-normal FLEET

%load in older data
load("WN2_FLEET.mat");
for i = 1:1
    PB_mean_velo = PB_FLEETResultsTable{:,1};
    PB_mean_velo_u = PB_FLEETResultsTable{:,2};
    PB_rms_velo = PB_FLEETResultsTable{:,3};
    PB_rms_velo_u = PB_FLEETResultsTable{:,4};
    PB_height = PB_FLEETResultsTable{:,10};
    PB_mod_velo = PB_FLEETResultsTable{:,12}; 
    PB_norm_SNR = PB_FLEETResultsTable{:,5};
    
    synth_mean_velo = synth_FLEETResultsTable{:,1};
    synth_mean_velo_u = synth_FLEETResultsTable{:,2};
    synth_rms_velo = synth_FLEETResultsTable{:,3};
    synth_rms_velo_u = synth_FLEETResultsTable{:,4};
    synth_height = synth_FLEETResultsTable{:,10};
    synth_mod_velo = synth_FLEETResultsTable{:,12};
    synth_norm_SNR = synth_FLEETResultsTable{:,5};
    
    oldPB_height = oldPB_FLEETResultsTable{:,10};
    oldPB_bin = oldPB_height<10;
    oldPB_height = oldPB_height(oldPB_bin);
    oldPB_mean_velo = oldPB_FLEETResultsTable{oldPB_bin,1};
    oldPB_mean_velo_u = oldPB_FLEETResultsTable{oldPB_bin,2};
    oldPB_rms_velo = oldPB_FLEETResultsTable{oldPB_bin,3};
    oldPB_rms_velo_u = oldPB_FLEETResultsTable{oldPB_bin,4};
    oldPB_norm_SNR = oldPB_FLEETResultsTable{:,5};
end

figure;
subplot(1,2,1);
plot(PB_mean_velo,PB_height,':m','Linewidth',2);
hold on;
    PB_upstream = 5491;
    for j = 1:length(PB_upstream)
        i = runs_list(uniqueruns==PB_upstream(j));
        plot(c_velocities{i},c_heights{i},'m','Linewidth',2);
    end
xlabel('$\bar{V}$','Interpreter','Latex')
ylabel('Height above surface [mm]');
legend(["Flat Plate PB","Inlet PB"])
set(gca,'FontSize', 18);
set(gca,'fontname','times')  % Set it to times
grid on;
ylim([0,15])
xlim([0,925])


subplot(1,2,2);
plot(synth_mean_velo,synth_height,':b','Linewidth',2);
hold on;
    SRA_upstream = [5492,5493];
    for j = 1:length(SRA_upstream)
        i = runs_list(uniqueruns==SRA_upstream(j));
        plot(c_velocities{i},c_heights{i},'b','Linewidth',2);
    end

xlabel('$\bar{V}$','Interpreter','Latex')
ylabel('Height above surface [mm]');
legend(["Flat Plate SRA","Inlet SRA"])
set(gca,'FontSize', 18);
set(gca,'fontname','times')  % Set it to times
grid on;
ylim([0,15])
xlim([0,925])

%% Testing reletive uncertainties
velo_edge = zeros(1,14);
velo_unc = zeros(1,14);
velo_RMS_unc = zeros(1,14);

for i = 1:14

velo_edge(i) = max(aw.u(:,locations_list(i)));
velo_unc(i) = median(unc_velo_COMBINED{i});
rms_unc = unc_rms_COMBINED{i};
velo_RMS_unc(i) = median(rms_unc(~isnan(rms_unc)));

end
% 
% mean(velo_RMS_unc)
% mean(velo_unc)
% 
% norm_unc_rms = mean(velo_RMS_unc./velo_edge)
% norm_uncmean = mean(velo_unc./velo_edge)

%% RMS plot
figure(50);
plot(DNS_RMS./max(DNS_V),DNS_h,'k','Linewidth',lnw)

        grid on;    
        xlabel('$V_{RMS}/V_e$','Interpreter','Latex')
        ylabel('Height above surface [mm]');
%         legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,.15]);
        ylim([0,im_height])
        hold on;

