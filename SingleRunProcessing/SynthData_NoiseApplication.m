function [noisey_synth_image] = SynthData_NoiseApplication(ideal_synth_image,...
    ROI_mean_SNR,ROI_gate1_location_bounds,ROI_gate2_location_bounds,...
    background_totalfit)
%This function applies noise to a synthetic image.

%% Identify rows in the first and second gate as a function of height
gatesmidpoint = mean([ROI_gate1_location_bounds(:,3),ROI_gate2_location_bounds(:,1)],2);
g1 = [round(ROI_gate1_location_bounds(:,1)-5),round(gatesmidpoint)];
g2 = [round(gatesmidpoint),round(ROI_gate2_location_bounds(:,3))+10];

%% Spottiness as f(y) and f(SNR)
noisey_synth_image = ones(size(ideal_synth_image));
ROI_mean_SNR(ROI_mean_SNR<1) = mean(ROI_mean_SNR);
g1(1:100,:) = repmat(g1(100,:),length(1:100),1);
g2(1:100,:) = repmat(g2(100,:),length(1:100),1);

for i = 1:length(ROI_mean_SNR)

    [spotiness_percent_g1] = GateSpotiness(1,0);

    spotiness_binary_g1 = randsrc(1,length(g1(i,1):g1(i,2)),[0,1;1-spotiness_percent_g1,spotiness_percent_g1]);
    noisey_synth_image(i,g1(i,1):g1(i,2)) = spotiness_binary_g1;

    [spotiness_percent_g2] = GateSpotiness(2,ROI_mean_SNR(i));

    spotiness_binary_g2 = randsrc(1,length(g2(i,1):g2(i,2)),[0,1;1-spotiness_percent_g2,spotiness_percent_g2]);
    noisey_synth_image(i,g2(i,1):g2(i,2)) = 1.6*spotiness_binary_g2;
end
    noisey_synth_image = noisey_synth_image.*ideal_synth_image;

%% Shot Noise
    scale_im = inv(50); %assuming an IRO amplification of the inverse of this number, like amplification of 500x. Larger amplification means worse SNR, less amplification means better SNR
    noisey_synth_image_scaled = noisey_synth_image.*scale_im;
    shotnoisey_synth_image_scaled = poissrnd(noisey_synth_image_scaled);
    shotnoiseonly_synth_image_scaled = noisey_synth_image_scaled-shotnoisey_synth_image_scaled;
    shotnoiseonly_synth_image = shotnoiseonly_synth_image_scaled./scale_im;
    noisey_synth_image_scaled = shotnoiseonly_synth_image+shotnoiseonly_synth_image;

%% Gaussian blurring
noisey_synth_image = imgaussfilt(noisey_synth_image,0.8);

%% Adding background and flare
noisey_synth_image = min(background_totalfit)+noisey_synth_image;
% 
%   figure;
%   imshow(noisey_synth_image,[1,round(max(noisey_synth_image(:)))]);

end