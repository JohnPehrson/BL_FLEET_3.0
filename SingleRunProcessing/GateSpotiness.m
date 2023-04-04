function [spotiness_percent] = GateSpotiness(gatenumber,SNR)
%This function will report the appropriate spotiness as a function of SNR.
%This is a function based on engineering judgement rather than some
%empirical formula. 
min_spotiness = 0.2;

max_snr = 35;
max_spotiness = 0.60;

    if gatenumber == 1
        spotiness_percent = 0.70;    
    else
        x_nondim = (SNR/max_snr)^2;
        spotiness_percent = x_nondim;
        spotiness_percent(spotiness_percent>max_spotiness) = max_spotiness;
    end
   spotiness_percent(spotiness_percent>1) = 1;
   spotiness_percent(spotiness_percent<=min_spotiness) = min_spotiness;

end