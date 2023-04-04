function [mean_velocity,rms_velocity,mean_velocity_r,mean_velocity_s,rms_velocity_r,rms_velocity_s,...
    mean_SNR,sufficient_counter] = Mean_RMS_Velo_Calculator(inst_velo,inst_velo_r,inst_velo_s,...
    SNR,numimages,filt_binary_master,emissionlocatingdata)
%This function calculates the mean and rms velocity

%% Initalizing Variables
allrows = 1:size(filt_binary_master,1);
mean_velocity = NaN(length(allrows),1);
rms_velocity = NaN(length(allrows),1);
mean_velocity_r = NaN(length(allrows),1);
mean_velocity_s = NaN(length(allrows),1);
rms_velocity_r = NaN(length(allrows),1);
rms_velocity_s = NaN(length(allrows),1);
mean_SNR = NaN(length(allrows),1);

%% Limiter to only gather data on rows with sufficient number of individual passed images
row_data_passed = sum(filt_binary_master,2)./numimages;
threshpass = 0.25;
test = row_data_passed>threshpass;
rows_ind = 1:length(row_data_passed);
test2 = bwareaopen(test, 50);
sufficient_counter = rows_ind(test2);

    %% Mean Velocity
    %Mean Velocity
    for i = sufficient_counter
        mean_velocity(i) = mean(inst_velo(i,filt_binary_master(i,:)));
    end
    
    %Mean Uncertainty
    for i = sufficient_counter
        N = sum(filt_binary_master(i,:));
        mean_velocity_r(i) = sqrt(sum((inst_velo_r(i,filt_binary_master(i,:)).^2)/(N^2)));
        mean_velocity_s(i) = sqrt(sum(inst_velo_s(i,filt_binary_master(i,:)).^2)/N);
    end
    
    % %Time-based velocity and centroid locations
    %     numimages = size(proc_velocity,1);
    %     timeaveraged_velocity = zeros(1,numimages);
    %     centroid_1 = zeros(1,numimages);
    %     centroid_2 = zeros(1,numimages);
    %     for i = 1:numimages
    %         if length(nonzeros(proc_velocity(i,:))) < 1
    %             timeaveraged_velocity(1,i) = 0;
    %             centroid_1(1,i) = 0;
    %             centroid_2(1,i) = 0;
    %         else
    %         timeaveraged_velocity(1,i) = mean(nonzeros(proc_velocity(i,:)));
    %         centroid_1(1,i) = mean(nonzeros(red_centroid_1(i,:)));
    %         centroid_2(1,i) = mean(nonzeros(red_centroid_2(i,:)));
    %         end
    %     end
    
    %% RMS Velocity
    %RMS
    deviation_velocity = inst_velo-mean_velocity;
    for i = sufficient_counter
        rms_velocity(i) = std(deviation_velocity(i,filt_binary_master(i,:)));
    end
    
    %RMS Uncertainty
    deviation_velocity_uncertainty_random = sqrt(inst_velo_r.^2+mean_velocity_r.^2);
    deviation_velocity_uncertainty_systematic = sqrt(inst_velo_s.^2+mean_velocity_s.^2);
    
    for i = sufficient_counter
        N = sum(filt_binary_master(i,:));
% Thesis version
    rms_velocity_r(i) = sqrt((N).*((1/(N*rms_velocity(i))).^2).*sum((deviation_velocity(i,filt_binary_master(i,:)).^2).*(deviation_velocity_uncertainty_random(i,filt_binary_master(i,:)).^2)));
    rms_velocity_s(i) = sqrt((N).*((1/(N*rms_velocity(i))).^2).*sum((deviation_velocity(i,filt_binary_master(i,:)).^2).*(deviation_velocity_uncertainty_systematic(i,filt_binary_master(i,:)).^2)));
    end
    
    %% R2 and signal and snr
    
    for i = sufficient_counter
        mean_SNR(i) = mean(SNR(i,filt_binary_master(i,:)));
    end

end

