clear all;close all;clc;

%% Time Average Comparison
%This script plots the time-averaged velocity and associated uncertainty
%from the fit to the time-averaged images. 

%% Overarching Variables
folder_path = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\ProcessedData\";
file_partial_name = "Time_Average_Fit_Run";
Run_Conditions_filepath = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/BLFLEETRunConditions.mat';   %stuff like gates and delays
Resolution_filepath = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/RefData.mat';                    %resolution
    load(Run_Conditions_filepath);
    load(Resolution_filepath);
runs_list = 1:14;
lens_standoff_dist = 300; %mm
height_focus = 3; %focusing around 3mm from the surface
decay_folderpath = "C:\Users\clark\Documents\GitHub\BL_FLEET_2.0\RawDataProcessing\ProcessedData\";
decay_filepaths = [ "ProcessedData_Synth_Run1_Images1500.mat";...
                    "ProcessedData_Synth_Run2_Images300.mat";...
                    "ProcessedData_Synth_Run3_Images300.mat"];

colorlist = ['r','b','g'];
PB_repeat   = [5484,5485,5486];
SRA_repeat  = [5479,5483,5495];
PB_BL       = [5491,5490,PB_repeat];
SRA_BL      = [5492,5493,5494,SRA_repeat];
PB_downstm  = [5487,5488,5489];

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

    %decay constant data sets
    c_taus          = cell(length(decay_filepaths),1);
    c_tau_errors    = cell(length(decay_filepaths),1);

%% Loading in Data
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
    c_wall_location_unc{i} = 1.5; %pixels
end

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
    res = res(1);
    run_centroids   = c_centroids{i};
    run_g1          = run_centroids(:,1);
    run_bounds_g1   = c_bounds{i,1};
    run_bounds_g2   = c_bounds{i,2};
    run_bounds_max  = run_bounds_g2(:,3)-run_bounds_g1(:,1);
    run_bounds_min  = run_bounds_g2(:,1)-run_bounds_g1(:,3);
    run_bounds_pix_max  = run_bounds_max.*res;
    run_bounds_pix_min  = run_bounds_min.*res;
    run_bounds_pix_max  = run_bounds_pix_max+run_centroids(:,1);
    run_bounds_pix_min  = run_bounds_pix_min+run_centroids(:,1);


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
%     
    %adjusting the actual data
    for i= 1:length(runs_list)
    
     %Put decay constant vector here instead of a constant------------------------------------------------
     decay_const_scalar = 450;
     decay_const_unc = [400,500];
     nondim_fit = polyval(p_tau_lin,decay_const_scalar);
     nondim_fit_unc = polyval(p_tau_lin,decay_const_unc);
     nondim_unc = abs(nondim_fit_unc(1)-nondim_fit_unc(2))/2;


     c_velocities{i} = (1+nondim_fit).*c_velocities{i};
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
            wall_loc_unc = c_wall_location_unc{i};
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

%% SNR Filtering for the max height
for i = 1:length(runs_list)
    
    %calculations to pick a maximum snr height
    snr_thresh = 3;
    no_filt_height = 7.5; %mm
    snr_binary = or((c_heights{i}<no_filt_height),(c_SNRs{i}'>snr_thresh)); %get data above some minimum height where rows have good snr
    above_surf_binary = (c_heights{i}>=0);
        h_mm = c_heights{i};
    height_snr_switch = h_mm(~snr_binary);
    height_lim = min(height_snr_switch);
    snr_binary2 = c_heights{i}<height_lim;
        both_binary = and(above_surf_binary,snr_binary2);



%     %plotting
%     figure(1);
%     subplot(1,3,1);
%     plot(c_velocities{i},c_heights{i},'b','Linewidth',2);
%     xlim([0,1000]);
%     ylim([0,max(c_heights{i})])
%     xlabel('Mean Velocity')
%     ylabel('Height above surface [mm]')
%     grid on;
%     set(gca,'FontSize', 15);
%     set(gca,'fontname','times')  % Set it to times
% 
%     subplot(1,3,2);
%     plot(c_SNRs{i},c_heights{i},'r','Linewidth',2);
%     xlabel('SNR')
%     grid on;
%     xlim([0,max(c_SNRs{i})]);
%     ylim([0,max(c_heights{i})])
%     hold off;
%     title(['SNR Filtering for Run ',num2str(uniqueruns(i))])
%     set(gca,'FontSize', 15);
%     set(gca,'fontname','times')  % Set it to times
%     xline(3,'k','Linewidth',3);
% 
%     subplot(1,3,3);
%     plot(snr_binary,c_heights{i},'r','Linewidth',2);
%     ylim([0,max(c_heights{i})])

    %filter data using SNR>threshold and a minimum height to keep
        f_velo          = c_velocities{i};
        f_velo_error    = c_velo_error_combined{i};
        f_height        = c_heights{i};
        f_velo          = f_velo(both_binary);
        f_velo_error    = f_velo_error(both_binary);
        f_height        = f_height(both_binary);

    c_velocities_plot{i}           = f_velo;
    c_velo_error_combined_plot{i}  = f_velo_error;
    c_heights_plot{i}              = f_height;
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
% legend('Mean Velocity','Mean Velocity 95% CI','Fitting Bounds');
set(gca,'FontSize', 15);
set(gca,'fontname','times')  % Set it to times
xlim([800,900]);
ylim([0,max(c_heights_plot{i})])
hold off;
end

%% Plotting relevant runs together
unc_linewidth = 1;

    %PB_repeatability
            labels_plot = strings(length(PB_repeat),1);
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
            end
        
            grid on;    
            title('Pizza Box Repeatability');
            xlabel('Velocity [m/s]');
            ylabel('Height above the surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([100,900]);
            ylim([0,max(c_heights_plot{i})])

     %SRA_repeatability
            labels_plot = strings(length(SRA_repeat),1);
            figure(4);
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
        
            grid on;    
            title('Synthetic Roughness Repeatability');
            xlabel('Velocity [m/s]');
            ylabel('Height above the surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([0,900]);
            ylim([0,max(c_heights_plot{i})])

    %SRA vs PB at 22C
            both_repeat = [PB_repeat,SRA_repeat];
            labels_plot = strings(length(both_repeat),1);
            figure(5);
            for j = 1:length(both_repeat)
                c = ceil(j/length(PB_repeat));
                i = runs_list(uniqueruns==both_repeat(j));
                plot(c_velocities_plot{i},c_heights_plot{i},colorlist(c),'Linewidth',2);
                hold on;
                labels_plot(j) = strcat("Run: ",num2str(both_repeat(j)));
            end
            for j = 1:length(both_repeat)
                i = runs_list(uniqueruns==both_repeat(j));
                c = ceil(j/length(PB_repeat));
                color_plot = [':',colorlist(c)];
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
            xlim([0,900]);
            ylim([0,max(c_heights_plot{i})])

    %PB_BL
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
        xlim([0,900]);
        ylim([0,max(c_heights_plot{i})])
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
        xlim([0,900]);
        ylim([0,max(c_heights_plot{i})])
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
            xlim([0,900]);
            ylim([0,max(c_heights_plot{i})])

    %SRA_BL
        j = 1;
        i = runs_list(uniqueruns==SRA_BL(j));
        color_plot = [':',colorlist(j)];
        downstream_loc_run = downstream_loc(i);

        figure(7);
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
        xlim([0,900]);
        ylim([0,max(c_heights_plot{i})])

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
        xlim([0,900]);
        ylim([0,max(c_heights_plot{i})])
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
        xlim([0,900]);
        ylim([0,max(c_heights_plot{i})])
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
            xlim([0,900]);
            ylim([0,max(c_heights_plot{i})])





