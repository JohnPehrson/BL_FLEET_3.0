clear all;close all;clc;

%% Time Average Comparison
%This script plots the time-averaged velocity and associated uncertainty
%from the fit to the time-averaged images. 

%% Overarching Variables
folder_path                 = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\ProcessedData\";
file_partial_name           = "Time_Average_Fit_Run";
Run_Conditions_filepath     = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/BLFLEETRunConditions.mat';   %stuff like gates and delays
Resolution_filepath         = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/RefData.mat';                    %resolution
ACE_On_Condition_filepath   = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/ACE_Data.mat';                    %resolution

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
cfd_filepath = "CFDDataForPlot.xlsx";


colorlist = ['r','b','g','m'];
PB_repeat   = [5484,5485,5486];
SRA_repeat  = [5479,5483,5495];
PB_BL       = [5491,5490,PB_repeat];
PB_BL_sp    = [5491,5490,5486];
SRA_BL      = [5492,5493,5494,SRA_repeat];
SRA_BL_sp   = [5493,5494,5495];
PB_downstm  = [5487,5488,5489]; %C, R, C-downstream

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

    %calculated/derived variables
    c_velocities                = cell(length(runs_list),1);
    c_heights_adjusted          = cell(length(runs_list),1);
    c_heights_unc               = cell(length(runs_list),1);
    c_random_velocity_error     = cell(length(runs_list),1);    %from fitting
    c_systematic_velocity_error = cell(length(runs_list),1);    %from fitting
    c_velo_error_combined       = cell(length(runs_list),1);
    c_heights_plot              = cell(length(runs_list),1);
    c_velocities_plot           = cell(length(runs_list),1);
    c_velo_error_combined_plot  = cell(length(runs_list),1);  
    c_velo_bounds               = cell(length(runs_list),1);  
    c_velo_error_wall_unc       = cell(length(runs_list),1);
    c_velo_error_magnification  = cell(length(runs_list),1);
    c_velo_error_span_point     = cell(length(runs_list),1);
    c_velo_error_decay          = cell(length(runs_list),1);
    c_SNR_binary_filt           = cell(length(runs_list),1);

    %decay constant data sets
    c_taus          = cell(length(decay_filepaths),1);
    c_tau_errors    = cell(length(decay_filepaths),1);

%% Loading in Data from the Sept 2022 Campaign
for i= 1:length(runs_list)
    load(fullfile(folder_path,ordered_filepaths(i)))
    c_centroids{i} = centroids;
    c_centroid_erros{i} = centroid_error;
    c_heights{i} = y_mm;
    c_SNRs{i} = snr;
    c_bounds(i,:) = {gb1,gb2};
    c_gates{i} = Gates(i,:);
    c_delays{i} = Delays(i,:);
    c_resolutions{i} = pixel_um_resolution(i,:);
    c_downstream_locs{i} = downstream_loc(i,:);
    c_downstream_locs_unc{i}= downstream_loc_unc(i,:);
    c_spanwise_locs{i} = spanwise_loc(i,:);
    c_spanwise_locs_unc{i} = spanwise_loc_unc(i,:);
    c_wall_location_unc{i} = zero_height_ref_unc; %pixels
    c_decay{i} = tau_fit;
end

%% Loading in Data from the (already processed) July 2022 Campaign
[Lam_FLEET,PB,SRA] = Load_July2022_Data(dist_LE_Trips,FLEET_July_2022_filepath);

%% Loading in Data from CFD that more closely matches the July 2022 Campaign
[Lam_CFD,DNS,RANS] = LoadCFD_Data_July2022(cfd_filepath,Run_Mean_Velo);
CFD_downstream_loc = 193; %mm
CFD_Reynolds = 5.22;  %Reynolds in Million/meter

%% Calculating Velocity
for i = 1:length(runs_list)
    run_cent        = c_centroids{i};
    run_cent_error  = c_centroid_erros{i};
    run_resolution  = c_resolutions{i};
    run_gates       = c_gates{i};
    run_delays      = c_delays{i};


    [velocity,random_velocity_error,systematic_velocity_error] = VelocityFinder(run_cent,...
        run_cent_error,run_resolution,run_gates,run_delays);

c_velocities{i}                 = velocity;
c_random_velocity_error{i}      = random_velocity_error;
c_systematic_velocity_error{i}  = systematic_velocity_error;
c_velo_error_combined{i} = sqrt(random_velocity_error.^2+systematic_velocity_error.^2);
end

%% Calculating min and max possible velocities if limited by fitting bounds
for i = 1:length(runs_list)
    res = c_resolutions{i};
    res = res(1)./1000;
    run_centroids   = c_centroids{i};
    run_g1          = run_centroids(:,1);
    run_bounds_g1   = c_bounds{i,1};
    run_bounds_g2   = c_bounds{i,2};
    run_bounds_max  = run_bounds_g2(:,3)-run_bounds_g1(:,1);
    run_bounds_min  = run_bounds_g2(:,1)-run_bounds_g1(:,3);
    run_bounds_pix_max  = run_bounds_max+run_centroids(:,1);
    run_bounds_pix_min  = run_bounds_min+run_centroids(:,1);

    run_cent_min    = [run_centroids(:,1),run_bounds_pix_min];
    run_cent_max    = [run_centroids(:,1),run_bounds_pix_max];
    run_cent_error  = [0,0]; %don't care
    run_resolution  = c_resolutions{i};
    run_gates       = c_gates{i};
    run_delays      = c_delays{i};


    [min_velo,~,~] = VelocityFinder(run_cent_min,run_cent_error,run_resolution,run_gates,run_delays);
    [max_velo,~,~] = VelocityFinder(run_cent_max,run_cent_error,run_resolution,run_gates,run_delays);

c_velo_bounds{i} = [min_velo,max_velo];
end


%% Finding the SNR Filtering for the max height
for i = 1:length(runs_list)
    %calculations to pick a maximum snr height
        snr_thresh = 3.15;
        no_filt_height = 7.5; %mm
        snr_binary = or((c_heights{i}<no_filt_height),(c_SNRs{i}'>snr_thresh)); %get data above some minimum height where rows have good snr
        above_surf_binary = (c_heights{i}>=0);
            h_mm = c_heights{i};
        height_snr_switch = h_mm(~snr_binary);
        height_lim = min(height_snr_switch);
        snr_binary2 = c_heights{i}<height_lim;
            both_binary = and(above_surf_binary,snr_binary2);


                %other run-specific filters
                if i==14
                    max_height = 11.25;
                    below_max_height = (c_heights{i}<=max_height);
                    both_binary = and(both_binary,below_max_height);
                end

        c_SNR_binary_filt{i} = both_binary;


        %plotting velo after filt
                f_velo          = c_velocities{i};
                f_height        = c_heights{i};
                f_velo          = f_velo(both_binary);
                f_height        = f_height(both_binary);
        


%             %plotting
%             figure(1);
%             subplot(1,3,1);
%             plot(c_velocities{i},c_heights{i},'b','Linewidth',2);
%             xlim([0,1000]);
%             ylim([0,max(c_heights{i})])
%             xlabel('Mean Velocity')
%             ylabel('Height above surface [mm]')
%             grid on;
%             set(gca,'FontSize', 15);
%             set(gca,'fontname','times')  % Set it to times
%         
%             subplot(1,3,2);
%             plot(c_SNRs{i},c_heights{i},'r','Linewidth',2);
%             xlabel('SNR')
%             grid on;
%             xlim([0,max(c_SNRs{i})]);
%             ylim([0,max(c_heights{i})])
%             hold off;
%             title(['SNR Filtering for Run ',num2str(uniqueruns(i))])
%             set(gca,'FontSize', 15);
%             set(gca,'fontname','times')  % Set it to times
%             xline(snr_thresh,'k','Linewidth',3);
%         
%             subplot(1,3,3);
%             plot(f_velo,f_height,'r','Linewidth',2);
%             ylim([0,max(c_heights{i})])
%             ylim([0,max(c_heights{i})])

end

%% Correcting for Emission Decay
    %get synthetic information from previous simulations
    for i = 1:length(decay_filepaths)
    load(fullfile(decay_folderpath,decay_filepaths(i)));

    binary_tau = synth_input_tau_fit<=440;

    c_taus{i} = synth_input_tau_fit(binary_tau);
    c_tau_errors{i} = nondim_velo_error(binary_tau);
    end

    %fit the data
    all_taus            = [c_taus{1};c_taus{2};c_taus{3}];
    all_nondim_errors   = [c_tau_errors{1};c_tau_errors{2};c_tau_errors{3}];

    plot_tau = linspace(min(all_taus),max(all_taus),100);
    p_tau_lin = polyfit(all_taus,all_nondim_errors,1);
    nondim_fit = polyval(p_tau_lin,plot_tau);

%     figure;
%     scatter(all_taus,all_nondim_errors);
%     hold on;
%     plot(plot_tau,nondim_fit);
    

    %adjusting the actual data
    for i= 1:length(runs_list)
    
        %fitting the decay constant of the data
            tau_run = c_decay{i};
            heights_run = c_heights{i};
            SNR_run = c_SNRs{i}';

        %use the snr binary filter to only get data with a good SNR
        %above the cutoff
            snr_filt = c_SNR_binary_filt{i};
            switchover_height_mm = 2;  %mm
            above_switchover_filt = heights_run>=switchover_height_mm;
            tau_filt = and(snr_filt,above_switchover_filt);
            
        %filter
            tau_run_f = tau_run(tau_filt);
            heigths_run_f = heights_run(tau_filt);

        %fit
            tau_fit_run = polyfit(heigths_run_f,tau_run_f,1);
            lin_fit_taus = polyval(tau_fit_run,heights_run);       
       %plotting fit   
%             figure;
%             scatter(tau_run,heights_run,'r');
%             hold on;
%             scatter(tau_run_f,heigths_run_f,'b');
%             plot(lin_fit_taus,heights_run,'k','Linewidth',2);
%             xlabel('Decay Constant [ns]');
%             ylabel('Height Above the Surface');
%             title(['Fitting tau as a line for Run ',num2str(uniqueruns(i))])
%                             
            %actually changing the real data
             decay_const_unc = [transpose(lin_fit_taus-50),transpose(lin_fit_taus+50)];
             nondim_fit = polyval(p_tau_lin,lin_fit_taus);
             nondim_fit_unc = polyval(p_tau_lin,decay_const_unc);
             nondim_unc = abs(nondim_fit_unc(:,1)-nondim_fit_unc(:,2))/2;
             nondim_fit = nondim_fit.*(0.75);  %%%%%%%%%%%%%%%  Get rid of this when I'm actually doing this for data I've fit using parameters from this dataset
        
                %metrics to manually check the amount of decay correction
                mean_decay_corr = mean(nondim_fit);
                min_max_decay_corr = [min(nondim_fit),max(nondim_fit)];
                disp(['Decay correction for Run ',num2str(uniqueruns(i)),' is '...
                    ,num2str(mean_decay_corr*100),'% with bounds of ',num2str(min_max_decay_corr(1).*100),'% to ',...
                    num2str(min_max_decay_corr(2).*100),'%'])

             c_velocities{i} = (1+nondim_fit').*c_velocities{i};
             c_velo_error_decay{i} = (nondim_unc).*c_velocities{i};
             c_velo_error_combined{i} = sqrt(c_velo_error_combined{i}.^2+c_velo_error_decay{i}.^2);
    end


%% Uncertainty Propogation (wall location, magnification, inclination angle, span angle, etc.)
for i= 1:length(runs_list)

    %     SpanwiseCameraAngle
    %     SpanwiseCameraAngle_unc %second column of SpanwiseCameraAngle

    %wall_loc_unc-------------------------------------------------------------
            heights = c_heights{i}';
            velo = c_velocities{i};
            wall_loc_unc = sqrt(c_wall_location_unc{i}.^2+1.^2);
            res = c_resolutions{i};
            [uncertainty_mean_velo_wall_loc,wall_loc_unc_mm] = Wall_Uncertainty_Propogator(res,...
                wall_loc_unc,velo,heights);
            c_heights_unc{i} = wall_loc_unc_mm.*ones(length(uncertainty_mean_velo_wall_loc),1);
            c_velo_error_wall_unc{i} = uncertainty_mean_velo_wall_loc;
            c_velo_error_combined{i} = sqrt(c_velo_error_combined{i}.^2+c_velo_error_wall_unc{i}.^2);

    %magnification------------------------------------------------------------
        L_hvec = heights-3; 
        inc_angle_max = inclination_angle_unc(i)+inclination_angle(i);
        angle_beam = 2; %likely beam angle, degrees
        d_vec = L_hvec.*sind(inc_angle_max+angle_beam);
        Mo_Md = 1-d_vec./lens_standoff_dist;
        mag_unc = abs(1-Mo_Md);
        c_velo_error_magnification{i} = c_velocities{i}.*mag_unc;
        c_velo_error_combined{i} = sqrt(c_velo_error_combined{i}.^2+c_velo_error_magnification{i}.^2);
        
    %Inclination angle
        c_heights_adjusted{i} = c_heights{i}'./cosd(inclination_angle(i));
        c_heights_unc{i} = sqrt(c_heights_unc{i}.^2+ abs(c_heights_adjusted{i}.*(sind(inclination_angle(i))./cosd(inclination_angle(i)))).*deg2rad(inclination_angle_unc(i)).^2);
        c_heights{i} = c_heights_adjusted{i}';

    %Span angle
        %span angle ideal is at 0 degrees, so mean and rms won't change, just
        %change systematic errors
        span_angle =  SpanwiseCameraAngle(i,1); %degree
        span_angle_uncertainty = SpanwiseCameraAngle(i,2); %degree
        %correct mean velocity
        c_velocities{i} = c_velocities{i}./cosd(span_angle);
        %account for uncertainty propogation
        c_velo_error_span_point{i} = sqrt((((c_velocities{i}/((cosd(span_angle).^2))).^2).*(deg2rad(span_angle_uncertainty).^2)));
        c_velo_error_combined{i} = sqrt(c_velo_error_combined{i}.^2+c_velo_error_span_point{i}.^2);
end

%% Smoothing the mean velocity and uncertainty for plotting
for i = 1:length(runs_list)
    r_velo      = c_velocities{i};
    r_velo_err  = c_velo_error_combined{i};
    r_height    = c_heights{i};

    r_velo(r_height<0) = 0;
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
%         plot(c_velocities{i}-c_velo_error_combined{i},c_heights{i},'k');
%         plot(c_velocities{i}+c_velo_error_combined{i},c_heights{i},'k');
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
    r_velo(replace_binary) = smooth_velo(replace_binary);
    r_velo_err(replace_binary) = smooth_velo_err(replace_binary);
    c_velocities{i}             = r_velo;
    c_velo_error_combined{i}    = r_velo_err;

end

%% Using the SNR filter to simplify the plotting, only plot where the SNR filter suggests
for i = 1:length(runs_list)

    %get filter
    snr_filt = c_SNR_binary_filt{i};

    %filter data using SNR>threshold and a minimum height to keep
        f_velo          = c_velocities{i};
        f_velo_error    = c_velo_error_combined{i};
        f_height        = c_heights{i};
        f_velo          = f_velo(snr_filt);
        f_velo_error    = f_velo_error(snr_filt);
        f_height        = f_height(snr_filt);

    c_velocities_plot{i}           = f_velo;
    c_velo_error_combined_plot{i}  = f_velo_error;
    c_heights_plot{i}              = f_height;


%   %plotting
%             figure(1);
%             subplot(1,3,1);
%             plot(c_velocities{i},c_heights{i},'b','Linewidth',2);
%             xlim([0,1000]);
%             ylim([0,max(c_heights{i})])
%             xlabel('Mean Velocity')
%             ylabel('Height above surface [mm]')
%             grid on;
%             set(gca,'FontSize', 15);
%             set(gca,'fontname','times')  % Set it to times
%         
%             subplot(1,3,2);
%             plot(c_SNRs{i},c_heights{i},'r','Linewidth',2);
%             xlabel('SNR')
%             grid on;
%             xlim([0,max(c_SNRs{i})]);
%             ylim([0,max(c_heights{i})])
%             hold off;
%             title(['SNR Filtering for Run ',num2str(uniqueruns(i))])
%             set(gca,'FontSize', 15);
%             set(gca,'fontname','times')  % Set it to times
%             xline(snr_thresh,'k','Linewidth',3);
%         
%             subplot(1,3,3);
%             plot(f_velo,f_height,'r','Linewidth',2);
%             ylim([0,max(c_heights{i})])
%             ylim([0,max(c_heights{i})])



end

%% Plotting Everything, in order
for i= 1:length(runs_list)

min_velo = c_velo_bounds{i};
min_velo = min_velo(:,1);
max_velo = c_velo_bounds{i};
max_velo = max_velo(:,2);

figure(2);
subplot(1,2,1);
plot(c_velocities{i},c_heights{i},'b','Linewidth',2);
xlim([0,1000]);
ylim([0,max(c_heights{i})])
xlabel('Mean Velocity')
ylabel('Height above surface [mm]')
grid on;
set(gca,'FontSize', 15);
set(gca,'fontname','times')  % Set it to times

subplot(1,2,2);
plot(c_velocities_plot{i},c_heights_plot{i},'b','Linewidth',2);
hold on;
plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},'k');

    plot(min_velo,c_heights{i},'r','Linewidth',2);
    plot(max_velo,c_heights{i},'r','Linewidth',2);

plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},'k');
grid on;    
title(['Mean Velocity for Run', num2str(uniqueruns(i))]);
xlabel('Velocity [m/s]');
ylabel('Height above the surface [mm]');
legend('Mean Velocity','Mean Velocity 95% CI','Fitting Bounds');
set(gca,'FontSize', 15);
set(gca,'fontname','times')  % Set it to times
xlim([0,900]);
ylim([0,max(c_heights_plot{i})])
hold off;
end

%% Plotting relevant runs together
unc_linewidth = 1;

    %PB_repeatability
            labels_plot = strings(length(PB_repeat),1);
            maxheight = 0;
            figure(3);
            for j = 1:length(PB_repeat)
                i = runs_list(uniqueruns==PB_repeat(j));
                plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(PB_repeat(j)));
            end
            for j = 1:length(PB_repeat)
                i = runs_list(uniqueruns==PB_repeat(j));
                color_plot = [':',colorlist(j)];
                plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                maxheight = max([maxheight,max(c_heights_plot{i})]);
            end
        
            grid on;    
            title('Pizza Box Repeatability');
            xlabel('Velocity [m/s]');
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
                plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(SRA_repeat(j)));
                maxheight = max([maxheight,max(c_heights_plot{i})]);
            end
            for j = 1:length(SRA_repeat)
                i = runs_list(uniqueruns==SRA_repeat(j));
                color_plot = [':',colorlist(j)];
                plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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

    %SRA vs PB at 22C
            both_repeat = [PB_repeat,SRA_repeat];
            labels_plot = ["PB Run","SRA Run"];
            proc_list = [1,4,2,3,5,6];
            maxheight = 0;
            figure(5);
            for j = proc_list
                c = ceil(j/length(PB_repeat));
                i = runs_list(uniqueruns==both_repeat(j));
                plot(c_velocities_plot{i},c_heights_plot{i},colorlist(c),'Linewidth',2);
                hold on;
%                 labels_plot(j) = strcat("Run: ",num2str(both_repeat(j)));
                maxheight = max([maxheight,max(c_heights_plot{i})]);
            end
            for j = 1:length(both_repeat)
                i = runs_list(uniqueruns==both_repeat(j));
                c = ceil(j/length(PB_repeat));
                color_plot = [':',colorlist(c)];
                plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
            end
        
            grid on;    
            title('Pizza Box vs Synthetic Roughness');
            xlabel('Velocity [m/s]');
            ylabel('Height above the surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([0,1000]);
            ylim([0,maxheight+1])

    %PB_BL seperate plots
            maxheight = 0;  
            for j= 1:length(PB_BL)
                i = runs_list(uniqueruns==PB_BL(j));
                maxheight = max([maxheight,max(c_heights_plot{i})]);
            end

        j = 1;
        i = runs_list(uniqueruns==PB_BL(j));
        color_plot = [':',colorlist(j)];
        downstream_loc_run = downstream_loc(i);

        figure(6);
        subplot(1,3,j);
        plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        grid on;    
        title(['PB ',num2str(downstream_loc_run),' mm']);
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        % legend('Mean Velocity','Mean Velocity 95% CI','Fitting Bounds');
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
        plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        grid on;    
        title(['PB ',num2str(downstream_loc_run),' mm']);
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        % legend('Mean Velocity','Mean Velocity 95% CI','Fitting Bounds');
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,1000]);
        ylim([0,maxheight+1])
        hold off;

        subplot(1,3,3)
        labels_plot = strings(length(PB_repeat),1);
            for j = 1:length(PB_repeat)
                i = runs_list(uniqueruns==PB_repeat(j));
                plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(PB_repeat(j)));
            end
            for j = 1:length(PB_repeat)
                i = runs_list(uniqueruns==PB_repeat(j));
                color_plot = [':',colorlist(j)];
                plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
            end
            
            downstream_loc_run = downstream_loc(i);
            grid on;    
            title(['PB ',num2str(downstream_loc_run),' mm']);
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
            plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
            hold on;
            labels_plot(j) = strcat(num2str(downstream_loc(i))," mm");
            maxheight = max([maxheight,max(c_heights_plot{i})]);
        end
        for j = 1:length(PB_BL_sp)
            i = runs_list(uniqueruns==PB_BL_sp(j));
            color_plot = [':',colorlist(j)];
            plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
            plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        end
    
        grid on;    
        title('PB Boundary Layer Evolution');
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,1000]);
        ylim([0,maxheight+1]);


    %SRA_BL
            maxheight = 0;  
            for j= 1:length(SRA_BL)
                i = runs_list(uniqueruns==SRA_BL(j));
                maxheight = max([maxheight,max(c_heights_plot{i})]);
            end

        j = 1;
        i = runs_list(uniqueruns==SRA_BL(j));
        color_plot = [':',colorlist(j)];
        downstream_loc_run = downstream_loc(i);

        figure(8);
        subplot(1,3,j);
        plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        grid on;    
        title(['SRA ',num2str(downstream_loc_run),' mm']);
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        % legend('Mean Velocity','Mean Velocity 95% CI','Fitting Bounds');
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,1000]);
        ylim([0,maxheight+1]);

        j = 2;
        i = runs_list(uniqueruns==SRA_BL(j));
        color_plot = [':',colorlist(j)];
        downstream_loc_run = downstream_loc(i);

        subplot(1,3,1);
        plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        grid on;    
        title(['SRA ',num2str(downstream_loc_run),' mm']);
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        % legend('Mean Velocity','Mean Velocity 95% CI','Fitting Bounds');
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
        plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
        hold on;
        plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        grid on;    
        title(['SRA ',num2str(downstream_loc_run),' mm']);
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        % legend('Mean Velocity','Mean Velocity 95% CI','Fitting Bounds');
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,1000]);
        ylim([0,maxheight+1]);
        hold off;

        subplot(1,3,3)
        labels_plot = strings(length(SRA_repeat),1);
            for j = 1:length(SRA_repeat)
                i = runs_list(uniqueruns==SRA_repeat(j));
                plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(SRA_repeat(j)));
            end
            for j = 1:length(SRA_repeat)
                i = runs_list(uniqueruns==SRA_repeat(j));
                color_plot = [':',colorlist(j)];
                plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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
                plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat(num2str(downstream_loc(i))," mm");
                maxheight = max([maxheight,max(c_heights_plot{i})]);
            end
            for j = 1:length(SRA_BL_sp)
                i = runs_list(uniqueruns==SRA_BL_sp(j));
                color_plot = [':',colorlist(j)];
                plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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

    %PB_downstm seperate subplots

        figure(10);
        maxheight = 0;
        labels_plot = strings(length(PB_downstm),1);
            for j = 1:length(PB_downstm)
                subplot(1,3,j);
                i = runs_list(uniqueruns==PB_downstm(j));
                plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(PB_downstm(j)));
                maxheight = max([maxheight,max(c_heights_plot{i})]);
            end
            for j = 1:length(PB_downstm)
                subplot(1,3,j);
                i = runs_list(uniqueruns==PB_downstm(j));
                color_plot = [':',colorlist(j)];
                plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
    
                downstream_loc_run = downstream_loc(i);
                spanwise_loc_run = spanwise_loc(i);
                grid on;    
                title(['PB ',num2str(downstream_loc_run),' mm downstream, ',num2str(spanwise_loc_run),' mm from center']);
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
            plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
            hold on;
            downstream_loc_run = downstream_loc(i);
            spanwise_loc_run = spanwise_loc(i);
            labels_plot(j) = strcat([num2str(downstream_loc_run),' mm downstream, ',num2str(spanwise_loc_run),' mm from center']);
        end
        for j = 1:length(PB_downstm)
            i = runs_list(uniqueruns==PB_downstm(j));
            color_plot = [':',colorlist(j)];
            plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
            plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
            maxheight = max([maxheight,max(c_heights_plot{i})]);
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

%% CFD Comparison
    %Pizza Box with  CFD, previous FLEET, and current FLEET all in one plot
    %current data
        labels_plot = strings(6,1);
        maxheight = 0;
        figure(12);
        for j = 1:length(PB_BL_sp)
            i = runs_list(uniqueruns==PB_BL_sp(j));
            plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
            hold on;
            labels_plot(j) = strcat(num2str(downstream_loc(i))," mm Downstream, FLEET, Walls On");
            maxheight = max([maxheight,max(c_heights_plot{i})]);
        end
%         for j = 1:length(PB_BL_sp)
%             i = runs_list(uniqueruns==PB_BL_sp(j));
%             color_plot = [':',colorlist(j)];
%             plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
%             plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
%         end

    %previous data
        plot(PB.Mean_Velo,PB.Heights,'m','Linewidth',2);
        labels_plot(4) = "193 mm Downstream, FLEET, Walls Off";
%         plot(PB.Mean_Velo-PB.Mean_Velo_Unc,PB.Heights,':m','Linewidth',unc_linewidth);
%         plot(PB.Mean_Velo+PB.Mean_Velo_Unc,PB.Heights,':m','Linewidth',unc_linewidth);

    %CFD
        plot(RANS.Velo,RANS.height,'k','Linewidth',2);
        plot(DNS.Velo,DNS.height,'--k','Linewidth',2);
        labels_plot(5) = "193 mm Downstream, RANS, Walls Off";
        labels_plot(6) = "193 mm Downstream, DNS, Periodic Domain";


        grid on;    
        title('PB Boundary Layer Evolution');
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,1000]);
        ylim([0,maxheight+1]);


%Synthetic Roughness with  CFD, previous FLEET, and current FLEET all in one plot
    %current data
        labels_plot = strings(6,1);
        maxheight = 0;
        figure(13);
        for j = 1:length(SRA_BL_sp)
            i = runs_list(uniqueruns==SRA_BL_sp(j));
            plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
            hold on;
            labels_plot(j) = strcat(num2str(downstream_loc(i))," mm Downstream, FLEET, Walls On");
            maxheight = max([maxheight,max(c_heights_plot{i})]);
        end
%         for j = 1:length(SRA_BL_sp)
%             i = runs_list(uniqueruns==PB_BL_sp(j));
%             color_plot = [':',colorlist(j)];
%             plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
%             plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
%         end

    %previous data
        plot(SRA.Mean_Velo,SRA.Heights,'m','Linewidth',2);
        labels_plot(4) = "193 mm Downstream, FLEET, Walls Off";
%         plot(PB.Mean_Velo-PB.Mean_Velo_Unc,PB.Heights,':m','Linewidth',unc_linewidth);
%         plot(PB.Mean_Velo+PB.Mean_Velo_Unc,PB.Heights,':m','Linewidth',unc_linewidth);

    %CFD
        plot(RANS.Velo,RANS.height,'k','Linewidth',2);
        labels_plot(5) = "193 mm Downstream, RANS, Walls Off";

        grid on;    
        title('SRA Boundary Layer Evolution');
        xlabel('Velocity [m/s]');
        ylabel('Height above the surface [mm]');
        legend(labels_plot(:));
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0,1000]);
        ylim([0,maxheight+1]);






