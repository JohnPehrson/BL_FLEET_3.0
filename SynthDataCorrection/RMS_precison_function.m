function [adj,rms_correction_unc] = RMS_precison_function(split_snr,synthrmsfit,...
    real_heights,halfCI,low_snr_mean,real,adj,low_snr_unc,coeffvals)
%This function provides the correction RMS as a function of SNR

    binary_snr = real.mean_SNR<split_snr;
    rms_correction = zeros(size(real.mean_SNR));
    rms_correction_unc = zeros(size(real.mean_SNR));
    for i = 1:length(real.mean_SNR)
        temp_snr = real.mean_SNR(i);
        height = real_heights(i);
        if binary_snr(i)   %low snr
            rms_correction(i) = synthrmsfit(temp_snr);
            up_b = halfCI+coeffvals;
            low_b = coeffvals-halfCI;
            up_cor = up_b(1)*(temp_snr.^up_b(2));
            low_cor = low_b(1)*(temp_snr.^low_b(2));
            rms_correction_unc(i) = abs(abs(low_cor)-abs(up_cor))/2;
        else %temp_snr<=split_snr %high snr
            rms_correction(i) = low_snr_mean;
            rms_correction_unc(i) = low_snr_unc;
        end
    end
    adj.velocity_rms = real.velocity_rms-rms_correction;
end