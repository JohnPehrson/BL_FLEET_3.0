clear all;close all;clc;

%% Time Average Comparison
%This script plots the time-averaged velocity and associated uncertainty
%from the fit to the time-averaged images. 

%% Overarching Variables
folder_path                 = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\ProcessedData\"; %"C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\ProcessedData\";
file_partial_name           = "Time_Average_Fit_RR1_Run";
Run_Conditions_filepath     = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/BLFLEETRunConditions.mat';   %stuff like gates and delays
Resolution_filepath         = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/RefData.mat';                    %resolution
ACE_On_Condition_filepath   = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/ACE_Data.mat';                    %resolution
near_wall_folder_path       = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\RunComparison\Near_Wall_Uncertainty_Data\";
near_wall_file_name         = "Time_Average_Fit_RR";
near_wall_folder_path_5479  = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\RunComparison\Near_Wall_Uncertainty_Data_5479\";

    load(Run_Conditions_filepath);
    load(Resolution_filepath);
    load(ACE_On_Condition_filepath);

runs_list = 1:14;
lens_standoff_dist = 300; %mm
N_ims = 1; %for time-averaged
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
all_runs    = [5479,5483:5495];

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
    c_heights_plot              = cell(length(runs_list),1);
    c_velocities_plot           = cell(length(runs_list),1);
    c_velo_bounds               = cell(length(runs_list),1);  
    c_SNR_binary_filt           = cell(length(runs_list),1);

    %Uncertainties
    unc_velo_centroid_fit       = cell(length(runs_list),1);    %from fitting
    unc_velo_resolution         = cell(length(runs_list),1);    %from fitting
    unc_height_resolution       = cell(length(runs_list),1);   
    unc_height_wall_location    = cell(length(runs_list),1); 
    unc_velo_emission_decay     = cell(length(runs_list),1); 
    unc_velo_magnification      = cell(length(runs_list),1); 
    unc_height_inclination      =cell(length(runs_list),1); 
    unc_velo_span_point         = cell(length(runs_list),1);
    unc_velo_nearwall           = cell(length(runs_list),1);
    unc_velo_RANDOM             = cell(length(runs_list),1);
    unc_velo_SYSTEMATIC         = cell(length(runs_list),1);
    unc_velo_COMBINED           = cell(length(runs_list),1);
    unc_height_SYSTEMATIC       = cell(length(runs_list),1);
    unc_velo_COMBINED_plot      = cell(length(runs_list),1); 
    unc_velo_heights            = cell(length(runs_list),1); 

    %Inner Variables
    inner_velo_fleet    = cell(length(runs_list),1); 
    inner_y_fleet       = cell(length(runs_list),1);
    inner_velo_cfd      = cell(length(runs_list),1);
    inner_y_cfd         = cell(length(runs_list),1); 

    %decay constant data sets
    c_taus          = cell(length(decay_filepaths),1);
    c_tau_errors    = cell(length(decay_filepaths),1);

    %save out data into an excel file
        %headers
        table_headers = {'Mean Velocity [m/s]','Uncertainty in Mean Velocity [m/s]','Height above the Test Article [mm]',...
        	            'Uncertainty in Height above the Test Article [mm]','Location Downstream the Leading Edge [mm]',...
                    	'Uncertainty in Downstream Location  [mm]','Spanwise Location [mm]','Uncertainty in Spanwise Location [mm]'};
        excel_filepath = "Mean_Velo_FLEET_Sept2022.xlsx";

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
    c_wall_location_unc{i} = zero_height_ref_unc; %pixels
    c_decay{i} = tau_fit;
end

%% Loading in Data from the (already processed) July 2022 Campaign
[Lam_FLEET,PB,SRA] = Load_July2022_Data(dist_LE_Trips,FLEET_July_2022_filepath);

%% Loading in Data from CFD that more closely matches the July 2022 Campaign
[Lam_CFD,DNS,RANS] = LoadCFD_Data_July2022(cfd_filepath,Run_Mean_Velo);
CFD_downstream_loc = 193; %mm
CFD_Reynolds = 5.22;  %Reynolds in Million/meter

%% Loading in Data for near-wall uncertainty 
[near_wall_fitobject] = NearWall_Uncertainty_Calculator(near_wall_folder_path,near_wall_file_name,Gates,Delays,pixel_um_resolution);
[near_wall_fitobject_5479] = NearWall_Uncertainty_Calculator(near_wall_folder_path_5479,near_wall_file_name,Gates,Delays,pixel_um_resolution);

%% Calculating Velocity
for i = 1:length(runs_list)
    run_cent        = c_centroids{i};
    run_cent_error  = c_centroid_erros{i};
    run_resolution  = c_resolutions{i};
    run_gates       = c_gates{i};
    run_delays      = c_delays{i};


    [velocity,centroid_fit_error,resolution_mean_velo_error] = VelocityFinder(run_cent,...
        run_cent_error,run_resolution,run_gates,run_delays);

c_velocities{i}                 = velocity;
unc_velo_centroid_fit{i}        = centroid_fit_error;
unc_velo_resolution{i}          = resolution_mean_velo_error;
% c_velo_error_combined{i}        = sqrt(centroid_fit_error.^2+resolution_mean_velo_error.^2);
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
        snr_thresh = 30;
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
        

            %plotting
            figure(1);
            subplot(1,3,1);
            plot(c_velocities{i},c_heights{i},'b','Linewidth',2);
            xlim([0,1000]);
            ylim([0,max(c_heights{i})])
            xlabel('Mean Velocity')
            ylabel('Height above surface [mm]')
            grid on;
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
        
            subplot(1,3,2);
            plot(c_SNRs{i},c_heights{i},'r','Linewidth',2);
            xlabel('SNR')
            grid on;
            xlim([0,max(c_SNRs{i})]);
            ylim([0,max(c_heights{i})])
            hold off;
            title(['SNR Filtering for Run ',num2str(uniqueruns(i))])
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xline(snr_thresh,'k','Linewidth',3);
        
            subplot(1,3,3);
            plot(f_velo,f_height,'r','Linewidth',2);
            ylim([0,max(c_heights{i})])
            ylim([0,max(c_heights{i})])

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
             unc_velo_emission_decay{i} = (nondim_unc).*c_velocities{i};
%              c_velo_error_combined{i} = sqrt(c_velo_error_combined{i}.^2+c_velo_error_decay{i}.^2);
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

    %magnification and beam angle------------------------------------------------------------
        L_hvec = heights-3; 
        inc_angle_max = inclination_angle_unc(i)+inclination_angle(i);
        angle_beam = 2; %likely beam angle, degrees
        d_vec = L_hvec.*sind(inc_angle_max+angle_beam);
        Mo_Md = 1-d_vec./lens_standoff_dist;
        mag_unc = abs(1-Mo_Md);
        unc_velo_magnification{i} = c_velocities{i}.*mag_unc;
        
    %Inclination angle
        c_heights{i} = c_heights{i}'./cosd(inclination_angle(i));
        unc_height_inclination{i} = c_heights{i}.*(sind(inclination_angle(i))./(cosd(inclination_angle(i).^2))).*deg2rad(inclination_angle_unc(i));

    %Span angle
        %change systematic errors
            span_angle =  SpanwiseCameraAngle(i,1); %degree
            span_angle_uncertainty = 2*SpanwiseCameraAngle(i,2); %degree
        %correct mean velocity
            c_velocities{i} = c_velocities{i}./cosd(span_angle);
        %account for uncertainty propogation
            unc_velo_span_point{i} = c_velocities{i}.*(sind(span_angle)./(cosd(span_angle.^2))).*deg2rad(span_angle_uncertainty);

    %Near-wall effects
    if i==1 %first run, different config
        unc_velo_nearwall{i} = feval(near_wall_fitobject_5479,c_heights{i});
    else
        unc_velo_nearwall{i} = feval(near_wall_fitobject,c_heights{i});
    end
end

%% Combine random and systematic uncertainties
for i = 1:length(runs_list)

    unc_velo_RANDOM{i} = (1/sqrt(N_ims)).*sqrt(unc_velo_centroid_fit{i}.^2);
    unc_velo_SYSTEMATIC{i} = sqrt(unc_velo_resolution{i}.^2+unc_velo_emission_decay{i}.^2+unc_velo_magnification{i}.^2 ...
                                +unc_velo_span_point{i}.^2+unc_velo_nearwall{i}.^2);
    unc_velo_COMBINED{i} = sqrt(unc_velo_RANDOM{i}.^2+unc_velo_SYSTEMATIC{i}.^2);
    unc_height_SYSTEMATIC{i} = sqrt(unc_height_resolution{i}.^2+unc_height_wall_location{i}.^2+unc_height_inclination{i}.^2);
end

%% Smoothing the mean velocity and uncertainty for plotting
for i = 1:length(runs_list)
    r_velo      = c_velocities{i};
    r_velo_err  = unc_velo_COMBINED{i};
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
    unc_velo_COMBINED{i}    = r_velo_err;

end

%% Propogating uncertainties in height into uncertainty in the mean velocity (mostly for plotting purposes)
unc_velo_COMBINED_table = unc_velo_COMBINED;
for i = 1:length(runs_list)
velo_unc    =  unc_velo_COMBINED{i};
height_unc  =  unc_height_SYSTEMATIC{i}; 
velo        =  c_velocities{i};
heights      =  c_heights{i};

        %move height up and down
            heights_dn = heights-height_unc;
            heights_up = heights+height_unc;
        %interpolate results if moved up or down
            velo_dn = interp1(heights_dn,velo,heights,'spline');
            velo_up = interp1(heights_up,velo,heights,'spline');
        %provide result
            uncertainty_mean_velo_wall_loc = (abs(velo_dn-velo)+abs(velo_up-velo))/2;
            unc_velo_heights{i} = uncertainty_mean_velo_wall_loc;
            unc_velo_COMBINED{i} = sqrt(unc_velo_COMBINED{i}.^2+unc_velo_heights{i}.^2);
    
            
%             %plot
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

%% Using the SNR filter to simplify the plotting and saving, only plot where the SNR filter suggests
for i = 1:length(runs_list)

    %get filter
    snr_filt = c_SNR_binary_filt{i};

    %filter data using SNR>threshold and a minimum height to keep
        f_velo              = c_velocities{i};
        f_velo_error        = unc_velo_COMBINED{i};
        f_velo_error_tbl    = unc_velo_COMBINED_table{i};
        f_height            = c_heights{i};
        f_height_unc        = unc_height_SYSTEMATIC{i};
        f_velo              = f_velo(snr_filt);
        f_velo_error        = f_velo_error(snr_filt);
        f_velo_error_tbl    = f_velo_error_tbl(snr_filt);
        f_height            = f_height(snr_filt);
        f_height_unc        = f_height_unc(snr_filt);

    c_velocities_plot{i}            = f_velo;
    unc_velo_COMBINED_plot{i}       = f_velo_error;
    unc_velo_COMBINED_table{i}      = f_velo_error_tbl;
    c_heights_plot{i}               = f_height;
    unc_height_SYSTEMATIC{i}        = f_height_unc;

end

% %% Plotting Everything, in order
% for i= 1:length(runs_list)
% 
% min_velo = c_velo_bounds{i};
% min_velo = min_velo(:,1);
% max_velo = c_velo_bounds{i};
% max_velo = max_velo(:,2);
% 
% figure(2);
% subplot(1,2,1);
% plot(c_velocities{i},c_heights{i},'b','Linewidth',2);
% xlim([0,1000]);
% ylim([0,max(c_heights{i})])
% xlabel('Mean Velocity')
% ylabel('Height above surface [mm]')
% grid on;
% set(gca,'FontSize', 15);
% set(gca,'fontname','times')  % Set it to times
% 
% subplot(1,2,2);
% plot(c_velocities_plot{i},c_heights_plot{i},'b','Linewidth',2);
% hold on;
% plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},'k');
% 
%     plot(min_velo,c_heights{i},'r','Linewidth',2);
%     plot(max_velo,c_heights{i},'r','Linewidth',2);
% 
% plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},'k');
% grid on;    
% title(['Mean Velocity for Run', num2str(uniqueruns(i))]);
% xlabel('Velocity [m/s]');
% ylabel('Height above the surface [mm]');
% legend('Mean Velocity','Mean Velocity 95% CI','Fitting Bounds');
% set(gca,'FontSize', 15);
% set(gca,'fontname','times')  % Set it to times
% xlim([0,900]);
% ylim([0,max(c_heights_plot{i})])
% hold off;
% end

%% Inner Variables
        Mach = Run_Mach(i,1);
        T_inf = 56; %K, static temp
        P_inf = 370; %Pa, static pressure
        U_inf = Run_Mean_Velo(i,1); %m/s
        use_runs = [PB_repeat,SRA_repeat];

    for i= 1:length(use_runs)
        j = runs_list(uniqueruns==use_runs(i));

        %FLEET
            y = c_heights_plot{j};
            u = c_velocities_plot{j};
            u_up = c_velocities_plot{j}+unc_velo_COMBINED_plot{j};
            u_dn = c_velocities_plot{j}-unc_velo_COMBINED_plot{j};
            [u_plus_c,u_vd_plus_c,y_plus_c] = InnerVariableCalculator(u,y,Mach,T_inf,P_inf,U_inf);
            [u_plus_c_up,u_vd_plus_c_up,~] = InnerVariableCalculator(u_up,y,Mach,T_inf,P_inf,U_inf);
            [u_plus_c_dn,u_vd_plus_c_dn,~] = InnerVariableCalculator(u_dn,y,Mach,T_inf,P_inf,U_inf);
                
        %CFD
            cfd_y = RANS.height;
            cfd_velo = RANS.Velo;
            U_inf = max(cfd_velo);
            [u_plus_c_cfd,u_vd_plus_c_cfd,y_plus_c_cfd] = InnerVariableCalculator(cfd_velo,cfd_y,Mach,T_inf,P_inf,U_inf);
    
        %Qualitative cutoff height 
        cutoff_height = 0.5; %mm
        [~,~,cutoff_inner_y] = InnerVariableCalculator(1000,cutoff_height,Mach,T_inf,P_inf,U_inf);

        %saving out data
            inner_velo_fleet{j} = [u_plus_c,u_plus_c_dn,u_plus_c_up];
            inner_y_fleet{j}    = y_plus_c;
            inner_velo_cfd{j}   = u_plus_c_cfd;
            inner_y_cfd{j}      = y_plus_c_cfd;
    end

%% Saving data into a table for CFD Validation
for i = 1:length(runs_list)

    velo            = round(c_velocities_plot{i},1);
    velo_unc        = round(unc_velo_COMBINED_table{i},1);
    heights         = round(c_heights_plot{i},3);
    heights_unc     = round(unc_height_SYSTEMATIC{i},3);
    down_loc        = downstream_loc(i).*ones(size(velo));
    down_loc_unc    = downstream_loc_unc(i).*ones(size(velo));
    span_loc        = spanwise_loc(i).*ones(size(velo));
    span_loc_unc    = spanwise_loc_unc(i).*ones(size(velo));

    T = table(velo,velo_unc,heights,heights_unc,down_loc,down_loc_unc,span_loc,span_loc_unc,'VariableNames',table_headers);
    sheet_name = ['Run ',num2str(uniqueruns(i))];
    writetable(T,excel_filepath,'Sheet',sheet_name,'Range','A1')

end

%% Plotting
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
                plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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
                plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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
                plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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
        plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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
        plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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
                plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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
            plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
            plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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
        plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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
        plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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
        plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
        plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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
                plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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
                plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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
                plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
    
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
            plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
            plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
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
        labels_plot = strings(5,1);
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
        


            %Plotting the most upstream SRA against the laminar profile
                %SRA stuff
                labels_plot = strings(length(SRA_BL(1:2)),1);
                maxheight = 0;
                figure(14);
                for j = 1:2
                    i = runs_list(uniqueruns==SRA_BL(j));
                    plot(c_velocities_plot{i},c_heights_plot{i},colorlist(j),'Linewidth',2);
                    hold on;
                    labels_plot(j) = strcat(num2str(c_downstream_locs{i}), ' mm Downstream, SRA FLEET, Walls On');
                    maxheight = max([maxheight,max(c_heights_plot{i})]);
                end
                %Laminar RANS
                    plot(Lam_CFD.Velo,Lam_CFD.height,'k','Linewidth',2);
                    labels_plot(3) = "193 mm Downstream, Laminar RANS, Walls Off";

                %previous data
                    plot(Lam_FLEET.Mean_Velo,Lam_FLEET.Heights,'m','Linewidth',2);
                    plot(Lam_FLEET.Mean_Velo-Lam_FLEET.Mean_Velo_Unc,Lam_FLEET.Heights,':m','Linewidth',unc_linewidth);
                    plot(Lam_FLEET.Mean_Velo+Lam_FLEET.Mean_Velo_Unc,Lam_FLEET.Heights,':m','Linewidth',unc_linewidth);
                    labels_plot(4) = "193 mm Downstream, Laminar FLEET, Walls Off";

                for j = 1:2
                    i = runs_list(uniqueruns==SRA_BL(j));
                    color_plot = [':',colorlist(j)];
                    plot(c_velocities_plot{i}-unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                    plot(c_velocities_plot{i}+unc_velo_COMBINED_plot{i},c_heights_plot{i},color_plot,'Linewidth',unc_linewidth);
                end
            
                grid on;    
                title('SRA vs Laminar');
                xlabel('Velocity [m/s]');
                ylabel('Height above the surface [mm]');
                legend(labels_plot(:));
                set(gca,'FontSize', 15);
                set(gca,'fontname','times')  % Set it to times
                xlim([0,1000]);
                ylim([0,maxheight+1])


%     %% Plotting all Runs
%             labels_plot = strings(length(all_runs),1);
%             maxheight = 0;
%             figure(16);
%             for j = 1:length(all_runs)
%                 i = runs_list(uniqueruns==all_runs(j));
%                 plot(c_velocities_plot{i},c_heights_plot{i},'Linewidth',2);
%                 hold on;
%                 labels_plot(j) = strcat("Run: ",num2str(all_runs(j)));
%             end
% %             for j = 1:length(all_runs)
% %                 i = runs_list(uniqueruns==all_runs(j));
% %                 plot(c_velocities_plot{i}-c_velo_error_combined_plot{i},c_heights_plot{i},'Linewidth',unc_linewidth);
% %                 plot(c_velocities_plot{i}+c_velo_error_combined_plot{i},c_heights_plot{i},'Linewidth',unc_linewidth);
% %                 maxheight = max([maxheight,max(c_heights_plot{i})]);
% %             end
% %         
%             grid on;    
%             title('All Runs');
%             xlabel('Velocity [m/s]');
%             ylabel('Height above the surface [mm]');
%             legend(labels_plot(:));
%             set(gca,'FontSize', 15);
%             set(gca,'fontname','times')  % Set it to times
%             xlim([0,1000]);
%             ylim([0,maxheight+1])

%% Inner Variable Plotting
ci_lw = 1;

%%______________________________Pizza Box_________________________%%

       %Ideal inner variable scaling
        y_ideal = logspace(1,2);
        xi = 0.41;
        C = 5;
        u_ideal = (1/xi).*log(y_ideal)+C;
            %viscous sublayer
        y_visc = logspace(-1,1);
        u_visc = y_visc;
                
        labels_plot = strings(7,1);
        figure(17);
        semilogx(y_ideal,u_ideal,'k','Linewidth',2);
        hold on;
        semilogx(y_visc,u_visc,'--k','Linewidth',2);
%         semilogx(y_plus_c_cfd,u_plus_c_cfd,'m','Linewidth',2);

        labels_plot(1:2) = ["Logarithmic Law of the Wall (0.41,5)";"Viscous Sublayer"];%;"RANS CFD"];
        
    for i= 1:length(PB_repeat)
        j = runs_list(uniqueruns==PB_repeat(i));

        %FLEET
            fleet_u = inner_velo_fleet{j};
            u_plus_c    = fleet_u(:,1);
            y_plus_c = inner_y_fleet{j};

        %plotting innver variable cfd
        semilogx(y_plus_c,u_plus_c,colorlist(i),'Linewidth',2);
        labels_plot(i+2) = strcat("FLEET - PB, Run ",num2str(uniqueruns(j)));
        
    end
    xline(cutoff_inner_y,':k','Linewidth',2)
    labels_plot(6) = "0.5 mm above surface";

  for i= 1:length(PB_repeat)
        j = runs_list(uniqueruns==PB_repeat(i));

        %FLEET
            fleet_u = inner_velo_fleet{j};
            u_plus_c    = fleet_u(:,1);
            u_plus_c_dn = fleet_u(:,2);
            u_plus_c_up = fleet_u(:,3);
            y_plus_c = inner_y_fleet{j};

        
        %plotting CI
                color_plot = [':',colorlist(i)];
        semilogx(y_plus_c,u_plus_c_up,color_plot,'Linewidth',ci_lw);
        semilogx(y_plus_c,u_plus_c_dn,color_plot,'Linewidth',ci_lw);
  end

    legend(labels_plot(:));
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times
    xlim([0.25,max(y_plus_c)+25])
    ylim([0,20])
    grid on;
    xlabel('y^+')
    ylabel('u^+_e_f_f')
    title('PB Inner Variables')

%%______________________________Synthetic Roughness_________________________%%

       %Ideal inner variable scaling
        y_ideal = logspace(1,2);
        xi = 0.41;
        C = 5;
        u_ideal = (1/xi).*log(y_ideal)+C;
            %viscous sublayer
        y_visc = logspace(-1,1);
        u_visc = y_visc;
                
        labels_plot = strings(7,1);
        figure(18);
        semilogx(y_ideal,u_ideal,'k','Linewidth',2);
        hold on;
        semilogx(y_visc,u_visc,'--k','Linewidth',2);
%         semilogx(y_plus_c_cfd,u_plus_c_cfd,'m','Linewidth',2);

        labels_plot(1:2) = ["Logarithmic Law of the Wall (0.41,5)";"Viscous Sublayer"];%;"RANS CFD"];
        
    for i= 1:length(SRA_repeat)
        j = runs_list(uniqueruns==SRA_repeat(i));

        %FLEET
            fleet_u = inner_velo_fleet{j};
            u_plus_c    = fleet_u(:,1);
            y_plus_c = inner_y_fleet{j};

        %plotting innver variable cfd
        semilogx(y_plus_c,u_plus_c,colorlist(i),'Linewidth',2);
        labels_plot(i+2) = strcat("FLEET - SRA, Run ",num2str(uniqueruns(j)));
    end
    xline(cutoff_inner_y,':k','Linewidth',2)
    labels_plot(6) = "0.5 mm above surface";

  for i= 1:length(SRA_repeat)
        j = runs_list(uniqueruns==SRA_repeat(i));

        %FLEET
            fleet_u = inner_velo_fleet{j};
            u_plus_c    = fleet_u(:,1);
            u_plus_c_dn = fleet_u(:,2);
            u_plus_c_up = fleet_u(:,3);
            y_plus_c = inner_y_fleet{j};

        
        %plotting CI
                color_plot = [':',colorlist(i)];
        semilogx(y_plus_c,u_plus_c_up,color_plot,'Linewidth',ci_lw);
        semilogx(y_plus_c,u_plus_c_dn,color_plot,'Linewidth',ci_lw);
  end

    legend(labels_plot(:));
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times
    xlim([0.25,max(y_plus_c)+25])
    ylim([0,20])
    grid on;
    xlabel('y^+')
    ylabel('u^+_e_f_f')
    title('SRA Inner Variables')
