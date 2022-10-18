%% Testing if I've got timing right;
load("Matfiles_preprocessing/timeaccurateunfiltereddata.mat");
[filt_centroids,filt_velocity,filt_velocity_s,filt_velocity_r,filt_R2,...
    filt_g2SNR,filt_signal,images_percentage_passed_filtering,filt_binary_master] = Filtering(red_centroids,...
    red_velocity,red_velocity_s, red_velocity_r,red_R2,red_g2SNR,red_signal,...
    gate1_location_bounds,gate2_location_bounds);


filt_binary_master(:,1:3000) = 0;
filt_binary_master(:,10000:end) = 0;

%% Mean
rows = transpose(size(filt_velocity,1):-1:1);
mean_velo = zeros(length(rows),1);
for i = 1:length(rows)
    row_velo = filt_velocity(i,filt_binary_master(i,:));
    mean_velo(i)= mean(row_velo);
end
figure;
plot(mean_velo,rows);