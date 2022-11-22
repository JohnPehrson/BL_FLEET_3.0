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

for i = 1:length(ROI_mean_SNR)

    [spotiness_percent_g1] = GateSpotiness(1,0);

    spotiness_binary_g1 = randsrc(1,length(g1(i,1):g1(i,2)),[0,1;1-spotiness_percent_g1,spotiness_percent_g1]);
    noisey_synth_image(i,g1(i,1):g1(i,2)) = spotiness_binary_g1;

    [spotiness_percent_g2] = GateSpotiness(2,ROI_mean_SNR(i));

    spotiness_binary_g2 = randsrc(1,length(g2(i,1):g2(i,2)),[0,1;1-spotiness_percent_g2,spotiness_percent_g2]);
    noisey_synth_image(i,g2(i,1):g2(i,2)) = spotiness_binary_g2;
end
    noisey_synth_image = noisey_synth_image.*ideal_synth_image;

%% Background noise
std_noise = mean(Amp_g2 ./ROI_mean_SNR(list_continuous));
r = randn(size(noisey_synth_image));
noise = 2.5.*(r.*std_noise+std_noise);
noise(noise<0) = 0;
noisey_synth_image = noisey_synth_image+noise;

low_data = prctile(noisey_synth_image(:),50);
bck_data = noisey_synth_image(noisey_synth_image(:)<low_data);
% std(bck_data(:))
% mean(bck_data(:))


%% Gaussian blurring
noisey_synth_image = imgaussfilt(noisey_synth_image,0.8);

% %% Plotting
% figure;
% image(noisey_synth_image);
% set(gca,'YDir','normal');
% colorbar;
% colormap(jet(round(max(noisey_synth_image(:)))));
% xlabel('Synth Streamwise Pixels')
% ylabel('Synth Vertical Pixels')
% title(['Single Synthetic Image'])
% grid on;


end