clear all;close all;clc;

%% Data visualizer, plotting intermediate data

runs = [1,2,3];
filepaths = ["ProcessedData/ProcessedData_Run1_Images929.mat";...
             "ProcessedData/ProcessedData_Run2_Images966.mat";...
             "ProcessedData/ProcessedData_Run3_Images748.mat"];
titles_runs = ["Laminar";"Pizza Box";"Synthetic Roughness"];
color_runs = ["b";"g";"r"];
color_runs_dotted = [":b";":g";":r"];

for i = runs(1):runs(end)
load(filepaths(runs(i)));

% mean_signal
% tau_fit
% mean_SNR
% velocity_mean
% velocity_rms
% velocity_mean_r
% velocity_mean_s
% velocity_rms_r
% velocity_rms_s
row_mm = velocimetry_geometricloc(:,5);
mean_uncertainty_combined = sqrt(velocity_mean_r.^2+velocity_mean_s.^2);
rms_uncertainty_combined =  sqrt(velocity_rms_r.^2+velocity_rms_s.^2);
max_height = max(row_mm);

figure(1)
hold on;
plot(velocity_mean,row_mm,color_runs(i),'Linewidth',2);
hold on;
% plot(velocity_rms,row_mm,'r','Linewidth',2);
% hold on;
plot(velocity_mean-mean_uncertainty_combined,row_mm,color_runs_dotted(i),'Linewidth',2);
hold on;
plot(velocity_mean+mean_uncertainty_combined,row_mm,color_runs_dotted(i),'Linewidth',2);
hold on;  
% plot(velocity_rms-rms_uncertainty_combined,row_mm,'k');
% hold on;
% plot(velocity_rms+rms_uncertainty_combined,row_mm,'k');
hold on;
grid on;    
ylim([0,max_height]);
title(['Mean Velocity for all Tripping Arrays']);
xlabel('Velocity [m/s]');
ylabel('Height above the surface [mm]');
% legend('Mean Velocity','RMS Velocity');
set(gca,'FontSize', 15);
set(gca,'fontname','times')  % Set it to times


figure(2)
hold on;
plot(velocity_rms,row_mm,color_runs(i),'Linewidth',2);
hold on;
% plot(velocity_rms-rms_uncertainty_combined,row_mm,color_runs_dotted(i),'Linewidth',2);
% hold on;
% plot(velocity_rms+rms_uncertainty_combined,row_mm,color_runs_dotted(i),'Linewidth',2);
grid on;    
ylim([0,max_height]);
title(['RMS Velocity for all Tripping Arrays']);
xlabel('Velocity [m/s]');
ylabel('Height above the surface [mm]');
% legend('Mean Velocity','RMS Velocity');
set(gca,'FontSize', 15);
set(gca,'fontname','times')  % Set it to times

figure(3)
hold on;
plot(mean_SNR,row_mm,color_runs(i),'Linewidth',2);
hold on;
grid on;    
ylim([0,max_height]);
title(['SNR for all Tripping Arrays']);
xlabel('SNR [-]');
ylabel('Height above the surface [mm]');
% legend('Mean Velocity','RMS Velocity');
set(gca,'FontSize', 15);
set(gca,'fontname','times')  % Set it to times

figure(4)
hold on;
plot(tau_fit,row_mm,color_runs(i),'Linewidth',2);
hold on;
grid on;    
ylim([0,max_height]);
title(['Decay Constants for all Tripping Arrays']);
xlabel('Decay Constant [ns]');
ylabel('Height above the surface [mm]');
% legend('Mean Velocity','RMS Velocity');
set(gca,'FontSize', 15);
set(gca,'fontname','times')  % Set it to times

    if i==1
        r1_heights = row_mm;
        r1_mean = velocity_mean;
        r1_mean_unc = mean_uncertainty_combined;
    elseif i==2
        r2_heights = row_mm;
        r2_mean = velocity_mean;
        r2_mean_unc = mean_uncertainty_combined;
    else
        r3_heights = row_mm;
        r3_mean = velocity_mean;
        r3_mean_unc = mean_uncertainty_combined;
    end


end

figure(2);
legend(titles_runs(1),titles_runs(2),titles_runs(3))

figure(3);
legend(titles_runs(1),titles_runs(2),titles_runs(3))

figure(4);
legend(titles_runs(1),titles_runs(2),titles_runs(3))


% figure;
% plot(mean_signal,row_mm,'b','Linewidth',2);

%% for legend
figure(5);
for i = 1:3
plot(mean_SNR,row_mm,color_runs(i),'Linewidth',2);
hold on;
set(gca,'FontSize', 15);
set(gca,'fontname','times')  % Set it to times
end
legend(titles_runs(1),titles_runs(2),titles_runs(3))





