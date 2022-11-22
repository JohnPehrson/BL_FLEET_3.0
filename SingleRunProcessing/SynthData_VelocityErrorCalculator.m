function [ideal_velocity_error,ideal_pixel_error,nondim_velo_error] = SynthData_VelocityErrorCalculator(gate2_trimmed,gate1_trimmed,gate2_image_ROI,...
            ideal_C2,ideal_C1,tau_fit,gate1_location_bounds,gate2_location_bounds,fitting_limits,...
            dt,scale,list_continuous,loc_midtime,heights_binary,velocity_mean)
%This function fits the two emission gates seperately and then compares
%the actual geometric centroid of the gates with the location the
%emission-fitting algorithm predicts is the centroid. This is then reported
%as a function of decay constant. 

%% Fit Gate 1 as f(y)
g1_fitcentroids = zeros(length(ideal_C2),1);
x = 1:size(gate2_trimmed,2);
rows = 1:size(gate2_trimmed,1);

    for i = rows
        fitting_limits(:,2) = gate1_location_bounds(i,:);
        fitting_limits(:,5) = gate2_location_bounds(i,:);
        x0 = fitting_limits(2,1:3);
        LB = fitting_limits(1,1:3);
        LB(3) = 2;
        UB = fitting_limits(3,1:3);
        y = gate1_trimmed(i,:);
        [temp_fitvariables] = SingleGaussFit(x0,LB,UB,x,y);
        g1_fitcentroids(i) = temp_fitvariables(2);
    end
    g1_diff = g1_fitcentroids-ideal_C1;


%% Fit Gate 2 as f(y)
 %fit the rows of data that aren't near the wall (minimal gate overlap)
g2_fitcentroids = zeros(length(ideal_C2),1);
x = 1:size(gate2_trimmed,2);
rows = 1:size(gate2_trimmed,1);

    for i = rows
        fitting_limits(:,2) = gate1_location_bounds(i,:);
        fitting_limits(:,5) = gate2_location_bounds(i,:);
        x0 = fitting_limits(2,4:6);
        LB = fitting_limits(1,4:6);
        LB(3) = 5;
        UB = fitting_limits(3,4:6);
        y = gate2_trimmed(i,:);
        [temp_fitvariables] = SingleGaussFit(x0,LB,UB,x,y);
        g2_fitcentroids(i) = temp_fitvariables(2);
    end
    g2_diff = g2_fitcentroids-ideal_C2;

%% Calculate velocity underprediction based on the centroid disparity
ideal_pixel_error = abs(g2_diff-g1_diff);
ideal_velocity_error = (ideal_pixel_error./(dt)).*scale(1);

verify_genmethod = ideal_C2-loc_midtime; %should be ~0

%% Plotting the location of centroids over the image
    figure;
    image(gate2_trimmed)
    colorbar;
    colormap(bone(3000));
    hold on;
    plot(gate2_location_bounds(:,1),rows,'r');
    hold on
    plot(g2_fitcentroids,rows,'b');
    hold on;
    plot(gate2_location_bounds(:,3),rows,'r');


%% Visualization
figure;
scatter(ideal_velocity_error(heights_binary),tau_fit(heights_binary));
title('Ideal Velocity Error vs Tau fit');
xlabel('Velocity Error [m/s]');
ylabel('Decay Constant');


%% Testing the velocity correction, to make sure the freestream is actually a line
tau_fit_eq = polyfit(tau_fit,ideal_velocity_error,1);
lin_fit_corr = polyval(tau_fit_eq,tau_fit);       
nondim_velo_error = 1+lin_fit_corr./max(prctile(velocity_mean,95));


       %plotting fit   
            figure;
            subplot(1,3,1);
            plot(velocity_mean,1:length(velocity_mean));
            set(gca, 'YDir','reverse')
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlabel('Mean Velocity');
            grid on;

            subplot(1,3,2);
            plot(tau_fit,1:length(velocity_mean));
            set(gca, 'YDir','reverse')
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlabel('Decay Constant');
            grid on;

            subplot(1,3,3);
            plot(nondim_velo_error.*velocity_mean,1:length(velocity_mean));
            set(gca, 'YDir','reverse')
            set(gca,'FontSize', 15);
            set(gca,'fontname','times')  % Set it to times
            xlabel('Corrected Mean Velocity');
            grid on;
            xline(880)


end