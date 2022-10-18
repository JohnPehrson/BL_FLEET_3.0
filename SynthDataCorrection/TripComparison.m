clear all;close all;clc;

%This function loads in the data for the three tripping cases. The goal is
%to compare the mean and rms velocity profiles. 

%% Filepaths
folderpath = 'C:\Users\clark\Documents\GitHub\BL_FLEET_2.0\SyntheticDataCorrection\CorrectedData';
runs = 1:4;  %Laminar, PB, Synthetic, old Pizza-Box
%% Load in CFD
cfdfilepath = "C:\Users\clark\Documents\GitHub\BL_FLEET_2.0\SyntheticDataCorrection\Case063-5-1.txt";
cfd = readtable(cfdfilepath);
cfd_y = cfd{:,3}.*1000;
cfd_velo = cfd{:,4};
frs_velo = [873.24]; %m/s

%% find file paths for the actual data
real_filepaths = strings(3,1);
a=dir(strcat(folderpath, "/*.mat"));
filenames = strings(length(a),1);
for i = 1:length(a)
filenames(i) = a(i).name;
end
for k = runs
    partial_input_filenames = strcat("CorrectedData_Run",num2str(runs(k)),"_Images");

    %find maximum number of images processed for real data
    real_filepaths(k) = filenames(contains(filenames,partial_input_filenames));
end

%% Load in tables
for k = runs
loadfilepath = fullfile(folderpath,real_filepaths(k));
load(loadfilepath)
    if k==1
    lam_FLEETResultsTable = FLEETResultsTable;
    elseif k==2
    PB_FLEETResultsTable = FLEETResultsTable;
    elseif k==3
    synth_FLEETResultsTable = FLEETResultsTable;
    else
    oldPB_FLEETResultsTable = FLEETResultsTable;
    end
end

%% Get data from the tables
lam_mean_velo = lam_FLEETResultsTable{:,1};
lam_mean_velo_u = lam_FLEETResultsTable{:,2};
lam_rms_velo = lam_FLEETResultsTable{:,3};
lam_rms_velo_u = lam_FLEETResultsTable{:,4};
lam_height = lam_FLEETResultsTable{:,10};
lam_mod_velo = lam_FLEETResultsTable{:,12};

PB_mean_velo = PB_FLEETResultsTable{:,1};
PB_mean_velo_u = PB_FLEETResultsTable{:,2};
PB_rms_velo = PB_FLEETResultsTable{:,3};
PB_rms_velo_u = PB_FLEETResultsTable{:,4};
PB_height = PB_FLEETResultsTable{:,10};
PB_mod_velo = PB_FLEETResultsTable{:,12}; 

synth_mean_velo = synth_FLEETResultsTable{:,1};
synth_mean_velo_u = synth_FLEETResultsTable{:,2};
synth_rms_velo = synth_FLEETResultsTable{:,3};
synth_rms_velo_u = synth_FLEETResultsTable{:,4};
synth_height = synth_FLEETResultsTable{:,10};
synth_mod_velo = synth_FLEETResultsTable{:,12};

oldPB_height = oldPB_FLEETResultsTable{:,10};
oldPB_bin = oldPB_height<10;
oldPB_height = oldPB_height(oldPB_bin);
oldPB_mean_velo = oldPB_FLEETResultsTable{oldPB_bin,1};
oldPB_mean_velo_u = oldPB_FLEETResultsTable{oldPB_bin,2};
oldPB_rms_velo = oldPB_FLEETResultsTable{oldPB_bin,3};
oldPB_rms_velo_u = oldPB_FLEETResultsTable{oldPB_bin,4};

%% Get data from CFD excel file
cfd_filepath = "CFDDataForPlot.xlsx";
DNS_V = readtable(cfd_filepath,'Range','E2:E202');
DNS_h = readtable(cfd_filepath,'Range','J2:J202');
Lam_V = readtable(cfd_filepath,'Range','O3:O203');
Lam_h = readtable(cfd_filepath,'Range','Q3:Q203');
RANS_V = readtable(cfd_filepath,'Range','V3:V203');
RANS_h = readtable(cfd_filepath,'Range','X3:X203');
DNS_V = DNS_V{:,:};
DNS_h = DNS_h{:,:};
Lam_V = Lam_V{:,:};
Lam_h = Lam_h{:,:};
RANS_V = RANS_V{:,:};
RANS_h = RANS_h{:,:};
scale_velo = frs_velo./max(Lam_V);
DNS_V = scale_velo.*DNS_V;
Lam_V = scale_velo.*Lam_V;
RANS_V = scale_velo.*RANS_V;

%% Plotting

%mean velocity laminar
figure;
hold on;
plot(lam_mean_velo,lam_height,'r','Linewidth',2);
% plot(lam_mean_velo+lam_mean_velo_u,lam_height,':r','Linewidth',2);
plot(lam_mean_velo-lam_mean_velo_u,lam_height,':r','Linewidth',2);

plot(Lam_V,Lam_h,'k','Linewidth',2);
plot(lam_mean_velo+lam_mean_velo_u,lam_height,':r','Linewidth',2);

ylim([0,16])
xlim([0,1000])
xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
ylabel('Height above the surface [mm]');
legend('FLEET','FLEET 95% CI','SA-RC-QCR CFD');
set(gca,'FontSize', 20);
set(gca,'fontname','times')  % Set it to times
title('FLEET vs RANS Laminar Case')
grid on;

%mean velocity 'turbulent'
figure;
% plot(RANS_V,RANS_h,'k','Linewidth',2);
hold on;
% plot(DNS_V,DNS_h,'r','Linewidth',2);

plot(lam_mean_velo,lam_height,'r','Linewidth',2);
plot(PB_mean_velo,PB_height,'g','Linewidth',2);
 plot(PB_mean_velo+PB_mean_velo_u,PB_height,':g','Linewidth',1.5);
% plot(PB_mean_velo-PB_mean_velo_u,PB_height,':g','Linewidth',2);


plot(synth_mean_velo,synth_height,'b','Linewidth',2);
plot(synth_mean_velo+synth_mean_velo_u,synth_height,':b','Linewidth',1.5);

plot(synth_mean_velo-synth_mean_velo_u,synth_height,':b','Linewidth',1.5);
plot(PB_mean_velo-PB_mean_velo_u,PB_height,':g','Linewidth',1.5);

ylim([0,16])
xlim([0,1000])
xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
ylabel('Height above the surface [mm]');
legend('Laminar','Pizza Box Array','Pizza Box 95% CI','Synthetic Roughness Array','Synthetic Roughness 95% CI');
set(gca,'FontSize', 20);
set(gca,'fontname','times')  % Set it to times
title('Mean Velocity')
grid on;

%RMS lam
figure;
plot(lam_rms_velo./frs_velo,lam_height,'r','Linewidth',2);
hold on;
plot(lam_rms_velo./frs_velo+lam_rms_velo_u./frs_velo,lam_height,':r','Linewidth',2);
hold on;
plot(lam_rms_velo./frs_velo-lam_rms_velo_u./frs_velo,lam_height,':r','Linewidth',2);

ylim([0,16])
xlim([0,0.1])
xlabel('$V_{RMS} / \bar{V}_{\infty}$ [m/s]','Interpreter','Latex')
ylabel('Height above the surface [mm]');
legend('RMS Velocity','95% CI');
set(gca,'FontSize', 20);
set(gca,'fontname','times')  % Set it to times
title('Laminar RMS Velocity')
grid on;

%RMS turb
figure;
plot(lam_rms_velo./frs_velo,lam_height,'r','Linewidth',2);
% plot(lam_mean_velo+lam_mean_velo_u,lam_height,':r','Linewidth',2);
% plot(lam_mean_velo-lam_mean_velo_u,lam_height,':r','Linewidth',2);
hold on;
plot(PB_rms_velo./frs_velo,PB_height,'g','Linewidth',2);
% plot(PB_mean_velo+PB_mean_velo_u,PB_height,':g','Linewidth',2);
% plot(PB_mean_velo-PB_mean_velo_u,PB_height,':g','Linewidth',2);

plot(synth_rms_velo./frs_velo,synth_height,'b','Linewidth',2);
% plot(synth_mean_velo+synth_mean_velo_u,synth_height,':b','Linewidth',2);
% plot(synth_mean_velo-synth_mean_velo_u,synth_height,':b','Linewidth',2);

ylim([0,16])
xlim([0,0.08])
xlabel('$V_{RMS} / \bar{V}_{\infty}$ [m/s]','Interpreter','Latex')
ylabel('Height above the surface [mm]');
legend('Laminar','Pizza Box','Synthetic Roughness');
set(gca,'FontSize', 20);
set(gca,'fontname','times')  % Set it to times
title('Root-Mean-Squared Velocity')
grid on;

%% Mini Plot
figure;
plot(RANS_V,RANS_h,'r','Linewidth',3);
ylim([0,16])
xlim([0,1000])
xlabel('$\bar{V}$ [m/s]','Interpreter','Latex')
ylabel('Height [mm]');
set(gca,'FontSize', 20);
set(gca,'fontname','times')  % Set it to times
grid on;

%% Inner variable scaling the FLEET measurements
    
    Mach = 5.73;
    T_inf = 56; %K, static temp
    P_inf = 400; %Pa, static pressure
    U_inf = 880; %m/s

    binary_synth_heights = synth_height<=13;
    binary_pb_heights = PB_height<=13;
    
    trim_synth_mean_velo = synth_mean_velo(binary_synth_heights);
    trim_synth_height = synth_height(binary_synth_heights);
    trim_PB_mean_velo = PB_mean_velo(binary_pb_heights);
    trim_PB_height = PB_height(binary_pb_heights);

    [u_plus_PB,u_vd_plus_PB,y_plus_PB,u_plus_c_PB,u_vd_plus_c_PB,...
        y_plus_c_PB] = InnerVariableCalculator(trim_PB_mean_velo,trim_PB_height,Mach,...
        T_inf,P_inf,U_inf);

    [u_plus_syn,u_vd_plus_syn,y_plus_syn,u_plus_c_syn,u_vd_plus_c_syn,...
        y_plus_c_syn] = InnerVariableCalculator(trim_synth_mean_velo,trim_synth_height,...
        Mach,T_inf,P_inf,U_inf);

%% Inner variable cfd, ideal, and plotting
    %load in all CFD data
    cfdfilepath = "C:\Users\clark\Documents\GitHub\BL_FLEET_2.0\SyntheticDataCorrection\CFDDataForPlot.xlsx";
    cfd = readtable(cfdfilepath);

    DNS_y_plus = cfd{:,2};
    DNS_u_vd_plus = cfd{:,6};
    RANS_y_plus = cfd{:,19};
    RANS_u_vd_plus = cfd{:,23};

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
        semilogx(y_ideal,u_ideal,'k','Linewidth',2);
        hold on;
        semilogx(y_visc,u_visc,'--k','Linewidth',2);
        hold on;
        semilogx(DNS_y_plus,DNS_u_vd_plus,'r','Linewidth',2);
        hold on;
        semilogx(RANS_y_plus,RANS_u_vd_plus,'g','Linewidth',2);
        hold on;
         semilogx(y_plus_c_PB,u_vd_plus_c_PB,'b','Linewidth',2);
        hold on;
        semilogx(y_plus_c_syn,u_vd_plus_c_syn,'m','Linewidth',2);

        legend('Logarithmic Law of the Wall (0.41,5)','Viscous Sublayer','DNS','RANS','Pizza Box','Synthetic Roughness')
        set(gca,'FontSize', 15);
        set(gca,'fontname','times')  % Set it to times
        xlim([0.4,500])
        ylim([0,20])
        grid on;
        xlabel('y^+')
        ylabel('u^+_e_f_f')
        title('FLEET vs CFD Inner Variables')









