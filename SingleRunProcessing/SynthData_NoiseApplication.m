function [noisey_synth_image] = SynthData_NoiseApplication(ideal_synth_image,...
    ROI_mean_SNR,ROI_gate1_location_bounds,ROI_gate2_location_bounds,...
    background_totalfit,Amp_g2,list_continuous)
%This function applies noise to a synthetic image.

%% Identify rows in the first and second gate as a function of height
gatesmidpoint = mean([ROI_gate1_location_bounds(:,3),ROI_gate2_location_bounds(:,1)],2);
g1 = [round(ROI_gate1_location_bounds(:,1)-5),round(gatesmidpoint)];
g2 = [round(gatesmidpoint),round(ROI_gate2_location_bounds(:,3))+10];

%% Spottiness as f(y) and f(SNR)
noisey_synth_image = ones(size(ideal_synth_image));
ROI_mean_SNR(ROI_mean_SNR<1) = min(ROI_mean_SNR(ROI_mean_SNR>0));
ROI_mean_SNR(end) = ROI_mean_SNR(end-1);
g1(1:100,:) = repmat(g1(100,:),length(1:100),1);
g2(1:100,:) = repmat(g2(100,:),length(1:100),1);

mean_spot_size_1 = 2;
std_spot_size_1 = 1;
mean_spot_size_2 = 3;
std_spot_size_2 = 3;

for i = 1:length(ROI_mean_SNR)

    [spotiness_percent_g1] = GateSpotiness(1,0);
        spotiness_binary_g1 = [];
        while length(spotiness_binary_g1)<(length(g1(i,1):g1(i,2))+20)
            subsequent_spot_size = round(normrnd( mean_spot_size_1, std_spot_size_1 ));
            spotiness_binary_g1 = [spotiness_binary_g1,(rand<spotiness_percent_g1)*ones(1,subsequent_spot_size)];
        end
        spotiness_binary_g1 = spotiness_binary_g1(10:(length(g1(i,1):g1(i,2))+9));
    noisey_synth_image(i,g1(i,1):g1(i,2)) = (1/spotiness_percent_g1).*spotiness_binary_g1;

    [spotiness_percent_g2] = GateSpotiness(2,ROI_mean_SNR(i));
        spotiness_binary_g2 = [];
        while length(spotiness_binary_g2)<(length(g2(i,1):g2(i,2))+20)
            subsequent_spot_size = round(normrnd( mean_spot_size_2, std_spot_size_2 ));
            spotiness_binary_g2 = [spotiness_binary_g2,(rand<spotiness_percent_g2)*ones(1,subsequent_spot_size)];
        end
        spotiness_binary_g2 = spotiness_binary_g2(10:(length(g2(i,1):g2(i,2))+9));
        dist_from_cent = ceil(-length(spotiness_binary_g2)/2:1:(length(spotiness_binary_g2)/2-1));
        dist_cent_mod = abs(dist_from_cent/max(dist_from_cent)).^2+1;
        noisey_synth_image(i,g2(i,1):g2(i,2)) = (1.75/spotiness_percent_g2).*spotiness_binary_g2.*dist_cent_mod;
    
end
    noisey_synth_image = noisey_synth_image.*ideal_synth_image;

%% Gaussian blurring
noisey_synth_image = imgaussfilt(noisey_synth_image,0.8);
    
%% Background noise
std_noise = mean(Amp_g2 ./ROI_mean_SNR(list_continuous));
r = randn(size(noisey_synth_image));
noise = 2.*(r.*std_noise+std_noise);
noise(noise<0) = 0;
noise_filt = imgaussfilt(noise,0.25);

noisey_synth_image = noisey_synth_image+noise_filt;
low_data = prctile(noisey_synth_image(:),50);
bck_data = noisey_synth_image(noisey_synth_image(:)<low_data);
% std(bck_data(:))
% mean(bck_data(:))
% 

%% Plotting
rows = 1:size(ideal_synth_image,1);
rows_bin = (rows<350) & (rows>84);
test = ideal_synth_image(rows_bin,:);
y = 1:size(test,1);
x = 1:size(test,2);


%real
figure;
image(ideal_synth_image(rows_bin,:));
set(gca,'YDir','reverse');
colormap(jet(3000));
% xlabel('Synth Streamwise Pixels')
% ylabel('Synth Vertical Pixels')
% title(['Ideal Synthetic Image'])
grid on;
axis off ; 
axis equal;
ylim([min(y),max(y)])
xlim([min(x),max(x)])

%synth
figure;
image(noisey_synth_image(rows_bin,:));
set(gca,'YDir','reverse');
colormap(jet(3000));
% xlabel('Synth Streamwise Pixels')
% ylabel('Synth Vertical Pixels')
% title(['Single Synthetic Image'])
grid on;
axis off ; 
axis equal;
ylim([min(y),max(y)])
xlim([min(x),max(x)])

savename = "Scitech_Pictures/Synth_Emissions.mat";
save(savename,'noisey_synth_image');



end