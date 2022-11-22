function [velocity_mean,heights_binary] = SynthSmoothVelo(velocity_mean,heights,run)
%This function smooths the mean velocity for synthetic data

    if run>=2 %aka either of the suspected turbulent cases
        %Measured data
        heights_binary = heights>0;
        heights_morezero = heights(heights_binary);
        
        %CFD Data
        T = readtable("C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\Case056NoTunnel.xlsx");
        T_height = T{:,1}.*1000;
        T_velo = T{:,2};
        
        %% Fitting the data up to the same height
        cfd_meanveloreplace = interp1(T_height,T_velo,heights_morezero);
        velocity_mean(heights_binary) = cfd_meanveloreplace;
        
        %plot
        figure;
        plot(T_height,T_velo);
        hold on;
        plot(heights,velocity_mean)
    else  %laminar
        velocity_mean;
        heights_binary = heights>0;
    end


end