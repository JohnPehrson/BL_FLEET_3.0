clear all;close all;clc;

%% Time Average Comparison
%This script plots the time-averaged velocity and associated uncertainty
%from the fit to the time-averaged images. 

%% Overarching Variables
folder_path                 = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\ProcessedData\"; %"C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\ProcessedData\";
file_partial_name           = "FullData_Run";
file_partial_name_synth     = "FullData_Synth_Run";
Run_Conditions_filepath     = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/BLFLEETRunConditions.mat';   %stuff like gates and delays
Resolution_filepath         = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/RefData.mat';                    %resolution
ACE_On_Condition_filepath   = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/ACE_Data.mat';                    %resolution
near_wall_folder_path       = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\RunComparison\Near_Wall_Uncertainty_Data\";
near_wall_file_name         = "Time_Average_Fit_RR";
near_wall_folder_path_5479  = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\RunComparison\Near_Wall_Uncertainty_Data_5479\";
inlet_cfd_filepath          = "Sept2022_RANS_CFD.xlsx";

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

colorlist = ['r','b','g','m'];
cfd_colorlist = ['k',':k'];
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
    synth_vert_offset        = [0,15,12,0,-5,0,0,0,0,0,0,-15,0,0];
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
    c_SNRs{i}                   = mean_SNR;
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
        ordered_filepaths = all_filenames;
%         ordered_filepaths = strings(length(uniqueruns),1);
%         for i = 1:length(ordered_filepaths)
%             filename_exp = strcat(file_partial_name_synth,num2str(runs_list(i)),"_");
%             which_run_binary = contains(all_filenames,filename_exp);
%             ordered_filepaths(i) = all_filenames(which_run_binary);
%         end

%j = 1:length(ordered_filepaths)

for j = 13
     
    load(fullfile(folder_path,ordered_filepaths(1)));

    s_velocities{j}             = velocity_mean;
    s_rms{j}                    = velocity_rms;
    s_rms_uncertainty{j}        = velocity_rms_s;
    s_SNRs{j}                   = mean_SNR;
    s_decay{j}                  = tau_fit;
    s_tau_input{j}              = synth_input_tau_fit;
    s_velo_corr{j}               = nondim_velo_error;
    s_heights{j}                = velocimetry_geometricloc(:,5);
end

%% Imprecision Correction using the RMS correction 
for j= 13 %1:length(runs_list)

    %measured data
    meas_rms =  c_rms{j};
    meas_snr =  c_SNRs{j};
    meas_rms_unc = unc_rms_resolution{j};

    %fit synthetic data
    snr = s_SNRs{j};
    rms = s_rms{j};
    unc = s_rms_uncertainty{j};
    heights = s_heights{j};
    heights_binary = heights>2.5;%mm
    snr = snr(heights_binary);
    rms = rms(heights_binary)+synth_vert_offset(j);
    unc = unc(heights_binary);
    synthrmsfit=fit(snr,rms,'power2'); 
    snr_stepper = linspace(min(snr),max(snr)+100,300);
       
    %calculate the corrected rms velocity
    rms_imprecision = synthrmsfit(meas_snr);
    for i = 1:length(rms_imprecision)
        if rms_imprecision(i) > meas_rms(i)
        rms_imprecision(i) = meas_rms(i);
        end
    end
    rms_real        = sqrt(meas_rms.^2-rms_imprecision.^2);
    rms_real_unc    = sqrt((meas_rms.^2).*(meas_rms_unc.^2)+(rms_imprecision.^2).*(mean(unc)))./rms_real;
        %lowest zero-value
        rms_row_stepper = 1:length(rms_real);
        zero_rows = rms_row_stepper(rms_real<5);
        use_below = max(zero_rows);
    rms_real(rms_row_stepper<=use_below) = NaN;
    rms_real_unc(rms_row_stepper<=use_below) = NaN;

    %not plotting rms velocity when uncertainty is especially large
    binary_not_rms = rms_real<rms_real_unc;
    rms_real(binary_not_rms) = NaN;

        %actually correct for rms velocity
        c_rms{j} = rms_real;
        unc_rms_synthetic{j} = rms_real_unc;

    

    figure;
    scatter(snr,rms,'b','filled');
    hold on;
    errorbar(snr,rms,unc,'b');
    scatter(meas_snr,meas_rms,'k','filled');
    plot(snr_stepper,synthrmsfit(snr_stepper),'r','Linewidth',2);
    scatter(meas_snr,rms_real,'k');
    errorbar(meas_snr,rms_real,rms_real_unc,'k');
    xlabel('SNR [-]');
    ylabel('URMS Velocity [m/s]');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    grid on;
    title('Imprecision estimate')
    ylim([0,max(rms)])
    legend('Synthetic','','Measured','Fit to Synthetic','Corrected RMS Velocity','');
    xlim([0,max(meas_snr)+10])

end

%% Correcting for Emission Decay
    %get synthetic information from previous simulations
    for i = 1:length(runs_list)
        %load data
            tau_fit             = s_decay{j};
            synth_input_tau_fit = s_tau_input{j};
            nondim_velo_error   = s_velo_corr{j};
        %fit data
            plot_tau = linspace(min(synth_input_tau_fit),max(synth_input_tau_fit),100);
            p_tau_lin = polyfit(synth_input_tau_fit,nondim_velo_error,1);
            nondim_fit_synth = polyval(p_tau_lin,plot_tau);

        %plot the fit of a line to the synthetic data
%             figure;
%             scatter(synth_input_tau_fit,nondim_velo_error);
%             hold on;
%             plot(plot_tau,nondim_fit_synth);
%             xlabel('Emission Decay Costant [ns]');
%             ylabel('Error in Mean Velocity');
%             set(gca,'FontSize', 15);
%             set(gca,'fontname','times')  % Set it to times
%             grid on;
%             title('Decay Correction')
%             legend('Synthetic Data','Fit to Synthetic Data');
 
        %fitting the decay constant of the data
            tau_run = c_decay{i};
            heights_run = c_heights{i};
            SNR_run = c_SNRs{i}';

        %use the snr binary filter to only get data with a good SNR
        %above the cutoff
            switchover_height_mm = 2;  %mm
            tau_filt = heights_run>=switchover_height_mm;
            
        %filter
            tau_run_f = tau_run(tau_filt);
            heigths_run_f = heights_run(tau_filt);

        %fit
            tau_fit_run = polyfit(heigths_run_f,tau_run_f,1);
            lin_fit_taus = polyval(tau_fit_run,heights_run);       
%        %plotting fit   
%             figure;
%             hold on;
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
             nondim_fit_measured    = polyval(p_tau_lin,lin_fit_taus);
             nondim_fit_unc         = polyval(p_tau_lin,decay_const_unc);
             nondim_unc = abs(nondim_fit_unc(:,1)-nondim_fit_unc(:,2))/2;
             nondim_fit_measured = ((nondim_fit_measured-1).*(0.75))+1;
             nondim_fit_measured(nondim_fit_measured<1) = 1;
        
                %metrics to manually check the amount of decay correction
                mean_decay_corr = mean(nondim_fit_measured)-1;
                min_max_decay_corr = [min(nondim_fit_measured)-1,max(nondim_fit_measured)-1];
                disp(['Decay correction for Run ',num2str(uniqueruns(i)),' is '...
                    ,num2str(mean_decay_corr*100),'% with bounds of ',num2str(min_max_decay_corr(1).*100),'% to ',...
                    num2str(min_max_decay_corr(2).*100),'%'])

             c_velocities{i} = (nondim_fit_measured').*c_velocities{i};
             unc_velo_emission_decay{i} = (nondim_unc).*c_velocities{i};
             c_rms{i} = (nondim_fit_measured').*c_rms{i};
             unc_rms_emission_decay{i} = (nondim_unc).*c_rms{i};
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
        unc_rms_magnification{i} = c_rms{i}.*mag_unc;
        
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
    unc_velo_RANDOM{i} = (1/sqrt(c_num_ims{i})).*sqrt(unc_velo_centroid_fit{i}.^2);
    unc_velo_SYSTEMATIC{i} = sqrt(unc_velo_resolution{i}.^2+unc_velo_emission_decay{i}.^2+unc_velo_magnification{i}.^2 ...
                                +unc_velo_span_point{i}.^2+unc_velo_nearwall{i}.^2);
    unc_velo_COMBINED{i} = sqrt(unc_velo_RANDOM{i}.^2+unc_velo_SYSTEMATIC{i}.^2);

    %rms velocity
    unc_rms_RANDOM{i} = (1/sqrt(c_num_ims{i})).*sqrt(unc_rms_centroid_fit{i}.^2);
    unc_rms_SYSTEMATIC{i} = sqrt(unc_rms_resolution{i}.^2+unc_rms_emission_decay{i}.^2+unc_rms_magnification{i}.^2 ...
                                +unc_rms_span_point{i}.^2+unc_rms_synthetic{i}.^2);
    unc_rms_COMBINED{i} = sqrt(unc_rms_RANDOM{i}.^2+unc_rms_SYSTEMATIC{i}.^2);

    %heights
    unc_height_SYSTEMATIC{i} = sqrt(unc_height_resolution{i}.^2+unc_height_wall_location{i}.^2+unc_height_inclination{i}.^2);
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

%% Plotting Mean Velocities
unc_linewidth = 1;

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
                labels_plot= [strcat("FLEET Run: ",num2str(uniqueruns(i))),"","","RANS CFD at matching Location"];
        grid on;    
        title(['PB ',num2str(downstream_loc_run),' mm']);
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
                labels_plot= [strcat("FLEET Run: ",num2str(uniqueruns(i))),"","","RANS CFD at matching Location"];
        grid on;    
        title(['PB ',num2str(downstream_loc_run),' mm']);
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
                labels_plot(j) = strcat("Run: ",num2str(PB_repeat(j)));
            end
                plot(aw.u(:,3),aw.h(:,3),cfd_colorlist(1),'Linewidth',2);
                labels_plot(4) = strcat("RANS CFD at matching Location");
            for j = 1:length(PB_repeat)
                i = runs_list(uniqueruns==PB_repeat(j));
                color_plot = [':',colorlist(j)];
                plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
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
            plot(c_velocities{i},c_heights{i},colorlist(j),'Linewidth',2);
            hold on;
            labels_plot(j) = strcat(num2str(downstream_loc(i))," mm");
            maxheight = max([maxheight,max(c_heights{i})]);
        end
        for j = 1:length(PB_BL_sp)
            i = runs_list(uniqueruns==PB_BL_sp(j));
            color_plot = [':',colorlist(j)];
            plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
            plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
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

%% Plotting RMS Velocities
unc_linewidth = 1;

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
            labels_plot(j) = strcat(num2str(downstream_loc(i))," mm");
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
lnw = 1;
cfd_lnw = 2;

figure(21); %men and rms velocity upstream of the shocks
v_inf = zeros(3,1);
subplot(2,3,1); %mean, 181 mm

        i = runs_list(uniqueruns==PB_BL(1));
        plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==SRA_BL(1));
        plot(c_velocities{i},c_heights{i},'r','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(2));
        plot(c_velocities{i},c_heights{i},'r','Linewidth',lnw);
        plot(aw.u(:,1),aw.h(:,1),cfd_colorlist(1),'Linewidth',cfd_lnw);
            v_inf(1) = max(aw.u(:,1));

    grid on;    
    title('Location 1');
    xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
    ylabel('Height above surface [mm]');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,900]);
    ylim([0,12]);

subplot(2,3,2); %mean, 320 mm

        i = runs_list(uniqueruns==PB_BL(2));
        plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==SRA_BL(3));
        plot(c_velocities{i},c_heights{i},'r','Linewidth',lnw);
        plot(aw.u(:,2),aw.h(:,2),cfd_colorlist(1),'Linewidth',cfd_lnw);
            v_inf(2) = max(aw.u(:,2));
        labels_plot= ["PB FLEET","SRA FLEET","RANS CFD"];

    legend(labels_plot)
    grid on;    
    title('Location 2');
    xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,900]);
    ylim([0,12]);

subplot(2,3,3); %mean, 370 mm

            i = runs_list(uniqueruns==PB_BL(3));
        plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==PB_BL(4));
        plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==PB_BL(5));
        plot(c_velocities{i},c_heights{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(4));
        plot(c_velocities{i},c_heights{i},'r','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(5));
        plot(c_velocities{i},c_heights{i},'r','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(6));
        plot(c_velocities{i},c_heights{i},'r','Linewidth',lnw);
        plot(aw.u(:,3),aw.h(:,3),cfd_colorlist(1),'Linewidth',cfd_lnw);
            v_inf(3) = max(aw.u(:,3));

    grid on;    
    title('Location 3');
    xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,900]);
    ylim([0,12]);

subplot(2,3,4); %RMS, 181 mm

            i = runs_list(uniqueruns==PB_BL(1));
        plot(c_rms{i}./v_inf(1),c_heights{i},'b','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==SRA_BL(1));
        plot(c_rms{i}./v_inf(1),c_heights{i},'r','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(2));
        plot(c_rms{i}./v_inf(1),c_heights{i},'r','Linewidth',lnw);

    grid on;    
    xlabel('$V_{RMS} / \bar{V}_{e}$','Interpreter','Latex')
    ylabel('Height above surface [mm]');
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,0.2]);
    ylim([0,6]);

subplot(2,3,5); %RMS, 320 mm

            i = runs_list(uniqueruns==PB_BL(2));
        plot(c_rms{i}./v_inf(2),c_heights{i},'b','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==SRA_BL(3));
        plot(c_rms{i}./v_inf(2),c_heights{i},'r','Linewidth',lnw);

    grid on;    
    xlabel('$V_{RMS} / \bar{V}_{e}$','Interpreter','Latex')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,0.2]);
    ylim([0,6]);

subplot(2,3,6); %RMS, 370 mm

            i = runs_list(uniqueruns==PB_BL(3));
        plot(c_rms{i}./v_inf(3),c_heights{i},'b','Linewidth',lnw);
        hold on;
            i = runs_list(uniqueruns==PB_BL(4));
        plot(c_rms{i}./v_inf(3),c_heights{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==PB_BL(5));
        plot(c_rms{i}./v_inf(3),c_heights{i},'b','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(4));
        plot(c_rms{i}./v_inf(3),c_heights{i},'r','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(5));
        plot(c_rms{i}./v_inf(3),c_heights{i},'r','Linewidth',lnw);
            i = runs_list(uniqueruns==SRA_BL(6));
        plot(c_rms{i}./v_inf(3),c_heights{i},'r','Linewidth',lnw);
            xline(0.5,':k');
    grid on;    
    xlabel('$V_{RMS} / \bar{V}_{e}$','Interpreter','Latex')
    set(gca,'FontSize', 15);
    set(gca,'fontname','times')  % Set it to times
    xlim([0,0.2]);
    ylim([0,6]);



    figure(22); %uncertainty upstream of the shocks
    
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
    xlim([0,0.15]);
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
                color_plot = [':',colorlist(j)];
                plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',unc_linewidth);
                    plot(aw.u(:,cfd_j(j)),aw.h(:,cfd_j(j)),cfd_colorlist(1),'Linewidth',2);
                downstream_loc_run = downstream_loc(i);
                spanwise_loc_run = spanwise_loc(i);
                grid on;    
                legend(strcat("Location ",num2str(j+3)'),'95% CI','','RANS CFD')
                xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
                ylabel('Height above surface [mm]');
                set(gca,'FontSize', 15);
                set(gca,'fontname','times')  % Set it to times
                xlim([0,900]);
                ylim([0,13]);
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
        ylim([0,6])


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
        for j = 1:length(both_repeat)
            i = runs_list(uniqueruns==both_repeat(j));
            if  ismember(both_repeat(j),PB_repeat) %pizza box
                c = 2;
            else
                c = 1;
            end
            color_plot = ['--',colorlist(c)];
            plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',1);
            plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',1);
        end
%         plot(aw.u(:,3),aw.h(:,3),cfd_colorlist(1),'Linewidth',2);

        clear labels_plot
        labels_plot = ["PB","","","PB 95% CI","","","","",""];
    
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

         both_repeat = [SRA_repeat(:)]; %3,1
        maxheight = 0;
        figure(27);
        subplot(1,2,2);
        for j = 1:length(both_repeat)
            if  ismember(both_repeat(j),PB_repeat) %pizza box
                c = 2;
            else
                c = 1;
            end
            i = runs_list(uniqueruns==both_repeat(j));
            plot(c_velocities{i},c_heights{i},colorlist(c),'Linewidth',1.5);
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
            plot(c_velocities{i}-unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',1);
            plot(c_velocities{i}+unc_velo_COMBINED{i},c_heights{i},color_plot,'Linewidth',1);
        end
%         plot(aw.u(:,3),aw.h(:,3),cfd_colorlist(1),'Linewidth',2);

        clear labels_plot
        labels_plot = ["SRA","","","SRA 95% CI","","","","",""];
    
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

%% Calculating t values for repeatability
points = 200;
h1 = zeros(points,3);
v1 = zeros(points,3);
s = zeros(points,3);
    for j = 1:length(PB_repeat)
        i = runs_list(uniqueruns==PB_repeat(j));
        v = c_velocities{i};
        h = c_heights{i};
        h1(:,j) = linspace(0,6.62839742571926,points);
        v1(:,j) = interp1(h,v,h1(:,j));
    end

    t_crit = tinv(0.90,3);
    h_mean = mean(h1,2);
    v_mean = mean(v1,2);

for k = 1:length(PB_repeat)
    for j = 1:length(h_mean)
        s(j,k) = (v_mean(j)-v1(j,k))/(t_crit/sqrt(3));
    end
end

points = 200;
h1 = zeros(points,3);
v1 = zeros(points,3);
s2 = zeros(points,3);
    for j = 1:length(SRA_repeat)
        i = runs_list(uniqueruns==SRA_repeat(j));
        v = c_velocities{i};
        h = c_heights{i};
        h1(:,j) = linspace(0,6.62839742571926,points);
        v1(:,j) = interp1(h,v,h1(:,j));
    end

    t_crit = tinv(0.90,3);
    h_mean = mean(h1,2);
    v_mean = mean(v1,2);

for k = 1:length(SRA_repeat)
    for j = 1:length(h_mean)
        s2(j,k) = (v_mean(j)-v1(j,k))/(t_crit/sqrt(3));
    end
end


    %Plot RMS velocity at L3 and downstream of the shocks on the same
    %plot
        mean_upstream_vinf = 862.457;
        figure(26);
                        i = runs_list(uniqueruns==PB_BL(3));
                        col = 'm';
                    plot(c_rms{i}./v_inf(3),c_heights{i},col,'Linewidth',lnw);
                    hold on;
                        i = runs_list(uniqueruns==PB_BL(4));
                    plot(c_rms{i}./v_inf(3),c_heights{i},col,'Linewidth',lnw);
                        i = runs_list(uniqueruns==PB_BL(5));
                    plot(c_rms{i}./v_inf(3),c_heights{i},col,'Linewidth',lnw);
                        i = runs_list(uniqueruns==SRA_BL(4));
                    plot(c_rms{i}./v_inf(3),c_heights{i},col,'Linewidth',lnw);
                        i = runs_list(uniqueruns==SRA_BL(5));
                    plot(c_rms{i}./v_inf(3),c_heights{i},col,'Linewidth',lnw);
                        i = runs_list(uniqueruns==SRA_BL(6));
                    plot(c_rms{i}./v_inf(3),c_heights{i},col,'Linewidth',lnw);
                        xline(0.5,':k');
                grid on;    
                xlabel('$V_{RMS} / \bar{V}_{e}$','Interpreter','Latex')
                set(gca,'FontSize', 15);
                set(gca,'fontname','times')  % Set it to times
                xlim([0,0.2]);
                ylim([0,6]);
    
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
        

            labels_plot = ["Location 3";"";"";"";"";"";"";labels_plot];
            grid on;    
            xlabel('$V_{RMS} / \bar{V}_{e}^*$','Interpreter','Latex')
            ylabel('Height above surface [mm]');
            legend(labels_plot(:));
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlim([0,0.2]);
            ylim([0,6])




