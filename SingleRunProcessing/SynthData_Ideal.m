function [ideal_synth_image,nondim_velo_error,synth_input_tau_fit,synth_input_velocity_mean,...
         ROI_mean_SNR,ROI_gate1_location_bounds,ROI_gate2_location_bounds,Amp_g2_mod,decay_error_eq...
         ] = SynthData_Ideal(rows,velocity_mean,heights,pixel_um_resolution,...
                                        Gates,Delays,doublegauss_fitvariables,...
                                        emissionlocatingdata,tau_fit,ROI,...
                                        gate1_location_bounds,gate2_location_bounds,fitting_limits,...
                                        list_continuous,mean_SNR,run)
%This function generates the ideal synthetic data intensity image
%replicating the real data except with known velocity information and with
%no velocity fluctuations

    %% Initialize image-wide information
    count_cols = length(ROI(3):ROI(4));
    count_rows = length(rows);
    %diffusion stuff
    p_w = [-10798366193.1642,252626.224750389,1.01596270692874];           
    heights_binary = heights>0;

    %% Calculate ideal displacement
    scale = pixel_um_resolution/(10^6);
    t = [ Delays(1)+Gates(1)/2 ,Delays(1)+Delays(2)+Gates(1)+Gates(2)/2];
    dt1 = (t(1))*10^(-9); %in seconds
    dx_real1 = (velocity_mean.*dt1)./(scale(1));
    dt2 = (t(2)-t(1))*10^(-9); %in seconds
    dx_real2 = (velocity_mean.*dt2)./(scale(1));

    %% Fitting information as f(y)
    amp_mod = 1;
    wid_mod = 0.9;
            Amp_g1 = amp_mod.*(doublegauss_fitvariables(:,1)/Gates(1)); 
            Width_g1 = wid_mod.*doublegauss_fitvariables(:,3);    %pixel width
            C1 = round(emissionlocatingdata(1))-max(dx_real1)+dx_real1;        %pixel location
                       
            Amp_g2 = amp_mod.*(doublegauss_fitvariables(:,4)/Gates(2));
            Width_g2 = doublegauss_fitvariables(:,6);    %pixel width
            ideal_C2 = C1+dx_real2;

    %% Synthetic Gate Generation
        %Models the gates used for FLEET emissions using amplitude decay,
        %diffusion, and fluid convection due to velocity. This will
        %generate both gates with one run

        [ideal_synth_image,gate2_trimmed,gate1_trimmed,gate2_image_ROI,...
         ROI_mean_SNR,ROI_gate1_location_bounds,ROI_gate2_location_bounds,...
         loc_midtime,Amp_g2_mod] = SynthGateGen(p_w,...
         Amp_g1,Width_g1,C1,Amp_g2,Width_g2,...
         count_cols,count_rows,emissionlocatingdata,...
         Delays,Gates,tau_fit,...
         velocity_mean,scale,ROI,mean_SNR,...
         gate1_location_bounds,gate2_location_bounds);
            
    %% Estimating the error in the mean velocity as a result of intensity decay
        [ideal_velocity_error,ideal_pixel_error,nondim_velo_error,decay_error_eq] = SynthData_VelocityErrorCalculator(gate2_trimmed,gate1_trimmed,gate2_image_ROI,...
            ideal_C2,C1,tau_fit,gate1_location_bounds,gate2_location_bounds,fitting_limits,dt2,scale,list_continuous,...
            loc_midtime,heights_binary,velocity_mean);
        
        synth_input_velocity_mean = velocity_mean;
        synth_input_tau_fit = tau_fit;
        
end