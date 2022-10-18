function [adj] = Syst_Precision(real,synth,adj,synth_binary,realheights_binary,...
                                folderpath,lam_real_filename,lam_real_stringlength,...
                                lam_synth_filename,lam_synth_stringlength,run)
%This function adjusts the systematic error in the rms velocity that arises
%due to precision errors

%% Fit the SNR near the wall as a line
% near_wall_binary = (real.velocimetry_geometricloc(:,5)>3)&(real.velocimetry_geometricloc(:,5)<10);
% wall_binary = real.velocimetry_geometricloc(:,5)<=3;
% wall_rows = real.velocimetry_geometricloc(wall_binary,5);
% near_wall_rows = real.velocimetry_geometricloc(near_wall_binary,5);
% near_wall_snr = real.mean_SNR(near_wall_binary);
% pf = polyfit(near_wall_rows,near_wall_snr,1);
% wall_snr = polyval(pf,wall_rows);
% real.mean_SNR(wall_binary) = wall_snr;
% 


%% Pre-plotting

figure;
scatter(real.mean_SNR,real.velocity_rms,'Linewidth',2);
hold on;

real_heights = real.velocimetry_geometricloc(:,5);

figure;
subplot(1,2,1);
plot(real.velocity_rms,real_heights,'r','Linewidth',2);
title('RMS Velocity')
xlabel('RMS Velocity [m/s]')
ylabel('Height above the surface [mm]')
set(gca,'FontSize', 15);
set(gca,'fontname','times')  % Set it to times
grid on;

subplot(1,2,2);
plot(real.mean_SNR,real_heights,'k','Linewidth',2);
title('SNR')
xlabel('SNR [-]')
set(gca,'FontSize', 15);
set(gca,'fontname','times')  % Set it to times
grid on;

%% Load in laminar data set for RMS subtraction
[lam_real,lam_synth,~] = RealSynthDataLoader(folderpath,...
lam_real_filename,lam_real_stringlength,lam_synth_filename,lam_synth_stringlength);
if run ~= 1
lam_real.mean_SNR = 1.2.*lam_real.mean_SNR;
end

%% Curve fitting the laminar data set in the portion where gate-overlap isn't a concern
height_bounds = [5,16];
heights_no_overlap  = (lam_real.velocimetry_geometricloc(:,5)>height_bounds(1))&(lam_real.velocimetry_geometricloc(:,5)<height_bounds(2));

x = lam_real.mean_SNR(heights_no_overlap);
y = lam_real.velocity_rms(heights_no_overlap);
synthrmsfit=fit(x,y,'power1'); 

snr_stepper = linspace(min(lam_real.mean_SNR),max(lam_real.mean_SNR),200);
figure;
scatter(lam_real.mean_SNR,lam_real.velocity_rms,'b','filled');
hold on;
plot(snr_stepper,synthrmsfit(snr_stepper),'r','Linewidth',2);
xlabel('SNR [-]');
ylabel('URMS Velocity [m/s]');
set(gca,'FontSize', 15);
set(gca,'fontname','times')  % Set it to times
grid on;
title('Precision estimate from Laminar Case')

%% Fitting the curve
offset_val = 8.8;
correction_factor = synthrmsfit(real.mean_SNR)-offset_val;
adj.velocity_rms = adj.velocity_rms-correction_factor;
%unceretainty
coeffvals= coeffvalues(synthrmsfit);
ci = confint(synthrmsfit,0.95);
halfCI = (ci(2,:)-ci(1,:))./2;
max_coeffs = coeffvals+[halfCI];
min_coeffs = coeffvals-[halfCI];
max_cf = max_coeffs(1).*(real.mean_SNR.^max_coeffs(2));
min_cf = min_coeffs(1).*(real.mean_SNR.^min_coeffs(2));
diff_cf = abs(max_cf-min_cf)./2;
adj.velocity_rms_s = sqrt(adj.velocity_rms_s.^2+diff_cf.^2);

%combined uncertainties
r_comb_unc = sqrt(real.velocity_rms_r.^2+real.velocity_rms_s.^2);
adj_comb_unc = sqrt(adj.velocity_rms_r.^2+adj.velocity_rms_s.^2);
r_heights = real.velocimetry_geometricloc(:,5);
adj_heights = adj.velocimetry_geometricloc(:,5);

%% Plotting
%precision error
figure;
scatter(real.mean_SNR,real.velocity_rms,'blue','Linewidth',2)
hold on;
scatter(adj.mean_SNR,adj.velocity_rms,'r','Linewidth',2)
title('Real vs Corrected RMS Velocity')
legend('Original Data','Corrected Data');
set(gca,'FontSize', 20);
set(gca,'fontname','times')  % Set it to times
xlabel('SNR [-]');
ylabel('RMS Velocity [m/s]');
grid on;

figure;
plot(real.velocity_rms,r_heights,'r','Linewidth',2)
hold on;
plot(adj.velocity_rms,adj_heights,'k','Linewidth',2)
hold on;
plot(real.velocity_rms-r_comb_unc,r_heights,':r','Linewidth',2)
hold on;
plot(real.velocity_rms+r_comb_unc,r_heights,':r','Linewidth',2)
hold on;
plot(adj.velocity_rms-adj_comb_unc,adj_heights,':k','Linewidth',2)
hold on;
plot(adj.velocity_rms+adj_comb_unc,adj_heights,':k','Linewidth',2)
title('Real vs Corrected RMS Velocity')
legend('Original Data','Corrected Data');
set(gca,'FontSize', 20);
set(gca,'fontname','times')  % Set it to times
xlabel('RMS Velocity [m/s]');
ylabel('Height above the surface [mm]');
grid on;


end