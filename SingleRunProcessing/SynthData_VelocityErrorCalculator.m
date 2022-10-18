function [ideal_velocity_error,ideal_pixel_error] = SynthData_VelocityErrorCalculator(gate2_trimmed,gate2_image_ROI,...
            ideal_C2,tau_fit,gate1_location_bounds,gate2_location_bounds,fitting_limits,...
            dt,scale,list_continuous,loc_midtime,heights_binary)
%This function fits only the second gate of FLEET emissions and compares
%the actual geometric centroid of the gate with the location the
%emission-fitting algorithm predicts is the centroid. This is then reported
%as a function of decay constant. 

%     %plotting
%     figure;
%     image(gate2_image_ROI)
%     colorbar;
%     maxplot = ceil(max(gate2_image_ROI(:)));
%     colormap(bone(maxplot));
%     hold on;
%     plot(gate1_location_bounds(:,1),list_continuous,'r');
%     hold on
%     plot(gate1_location_bounds(:,2),list_continuous,'b');
%     hold on;
%     plot(gate1_location_bounds(:,3),list_continuous,'r');
%     hold on;
%     plot(gate2_location_bounds(:,1),list_continuous,'r');
%     hold on
%     plot(gate2_location_bounds(:,2),list_continuous,'b');
%     hold on;
%     plot(gate2_location_bounds(:,3),list_continuous,'r');


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
% 
% 
%         if i>50
%         fit_gauss = @(p) p(1)*exp(-(x-p(2)).^2./(2*p(3)^2));
%         figure;
%         plot(x,y,'k');
%         hold on;
%         plot(x,fit_gauss(temp_fitvariables),'--r');
%         hold on;
%         plot([ideal_C2(i),ideal_C2(i)],[0,2000],'k')
%         hold on;
%         plot([g2_fitcentroids(i),g2_fitcentroids(i)],[0,2000],'r')
%         legend('Ideal Gate Intensity','Fit','Middle of Emissions','Fit Centroid')
%         disp('test');
%         end


    end


%% Calculate velocity underprediction based on the centroid disparity
ideal_pixel_error = abs(ideal_C2-g2_fitcentroids);
ideal_velocity_error = (ideal_pixel_error./(dt)).*scale(1);

verify_genmethod = ideal_C2-loc_midtime; %should be ~0
%loc_midtime

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


% %% Visualization
% figure;
% scatter(ideal_velocity_error(heights_binary),tau_fit(heights_binary));
% title('Ideal Velocity Error vs Tau fit');
% xlabel('Velocity Error [m/s]');
% ylabel('Decay Constant');
end