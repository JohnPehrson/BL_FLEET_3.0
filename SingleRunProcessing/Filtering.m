function [out_centroids,out_velocity,out_velocity_s,out_velocity_r,...
    out_g2SNR,images_percentage_passed_filtering,filt_binary_master_passed] = Filtering(red_centroids,...
    red_velocity,red_velocity_s,red_velocity_r,red_g2SNR,gate1_location_bounds,gate2_location_bounds,...
    synth_switch,row_mm)
%This sub-function filters various parameters from the fitting to make sure
%only data with sufficient quality is used for final data calculations

%% Initializing variables
g1_centroids = red_centroids(:,1:2:end);
g2_centroids = red_centroids(:,2:2:end);
out_centroids = zeros(size(red_centroids));

%% Setting thresholds for the various parameters and creating a binary map
%centroids not near end bounds
nearbounds = 1;
filt_g1_bounds_input = and(abs(gate1_location_bounds(:,1)-g1_centroids)>nearbounds,abs(gate1_location_bounds(:,3)-g1_centroids)>nearbounds);
filt_g2_bounds_input = and(abs(gate2_location_bounds(:,1)-g2_centroids)>nearbounds,abs(gate2_location_bounds(:,3)-g2_centroids)>nearbounds);
filt_g1_bounds = filt_g1_bounds_input;
filt_g2_bounds = filt_g2_bounds_input;

numims = size(filt_g1_bounds_input,2);
rows = size(filt_g1_bounds_input,1);
row_dist_filt = 1; %how close filtered row must be to also filter nearby rows
for i = 1:numims
    for j = 2:(rows-1)
        filt_g1_bounds(j,i) = filt_g1_bounds_input(j-row_dist_filt,i)&filt_g1_bounds_input(j,i)&filt_g1_bounds_input(j+row_dist_filt,i); %%if nearby row is filtered, also filter this row
        filt_g2_bounds(j,i) = filt_g2_bounds_input(j-row_dist_filt,i)&filt_g2_bounds_input(j,i)&filt_g2_bounds_input(j+row_dist_filt,i); %%if nearby row is filtered, also filter this row
    end
end

%filtering rows with abnormally low SNR (rowwise filtering)
snr_thresh = 1; %# of std away from mean
mean_row_snr = mean(red_g2SNR,2);
std_row_snr = std(red_g2SNR')';
lower_bound_snr = mean_row_snr-snr_thresh.*std_row_snr;
lower_bound_snr(lower_bound_snr<0) = 0;
filt_SNR_rowwise = red_g2SNR>lower_bound_snr;

%sufficient SNR based on image average
snr_thresh_min = 6; %# minimum permitted average SNR per row
mean_row_signal = mean(red_g2SNR,2);
filt_SNR = mean_row_signal>snr_thresh_min;

%data with sufficiently small random error/uncertainty
r_velo_thresh = 1.2; %# of std away from mean
mean_r_velo = mean(red_velocity_r,2);
std_row_r_velo = std(red_velocity_r')';
lower_bound_r_velo = mean_r_velo-r_velo_thresh.*std_row_r_velo;
lower_bound_r_velo(lower_bound_r_velo<0) = 0;
lower_bound_r_velo(lower_bound_r_velo>250) = 250;
filt_r_velo_rowwise = red_velocity_r>lower_bound_r_velo;

% filtering using synthetic data, limit upper bound
if synth_switch
    filt_synth = ones(size(filt_r_velo_rowwise))>0;
    temp_filt_rows = row_mm<14;
    filt_synth = filt_synth.*temp_filt_rows;
else
    filt_synth = ones(size(filt_r_velo_rowwise))>0;
end

%combine all filters into the master filter
filt_binary_master = ~(filt_g2_bounds & filt_SNR & filt_r_velo_rowwise & filt_SNR_rowwise & filt_synth); %& filt_r2
filt_binary_master_passed = ~filt_binary_master;
images_percentage_passed_filtering = 1-mean(filt_binary_master,2);

figure;
plot(images_percentage_passed_filtering,1:length(images_percentage_passed_filtering));
hold on;
plot(mean(filt_g2_bounds,2),1:length(images_percentage_passed_filtering));
hold on;
plot(mean(filt_SNR,2),1:length(images_percentage_passed_filtering));
hold on;
plot(mean(filt_SNR_rowwise,2),1:length(images_percentage_passed_filtering));
hold on;
plot(mean(filt_r_velo_rowwise,2),1:length(images_percentage_passed_filtering));
hold on;
legend('Total Filtering by row','G2 Bounds Filtering','SNR Filtering Threshold','SNR Filtering Row-wise','Velocity Uncertainty (Random) Filtering');
set(gca, 'YDir','reverse')

%% Filtering data in output variable format
g1_centroids(filt_binary_master) = NaN;
g2_centroids(filt_binary_master) = NaN;

red_velocity(filt_binary_master) = NaN;
red_velocity_s(filt_binary_master) = NaN;
red_velocity_r(filt_binary_master) = NaN;
red_g2SNR(filt_binary_master) = NaN;


out_centroids(:,1:2:end) = g1_centroids;
out_centroids(:,2:2:end) = g2_centroids;
out_velocity = red_velocity;
out_velocity_s = red_velocity_s;
out_velocity_r = red_velocity_r;
out_g2SNR = red_g2SNR;
end