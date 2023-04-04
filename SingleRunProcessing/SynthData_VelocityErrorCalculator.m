function [ideal_velocity_error,ideal_pixel_error,nondim_velo_error,decay_error_eq] = SynthData_VelocityErrorCalculator(gate2_trimmed,gate1_trimmed,gate2_image_ROI,...
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
fitting_limits(:,6) = [4;6;8];

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

        if i==200
            x_dense = linspace(min(x),max(x),1000);
            gaussdist = @(p) p(1)*exp(-(x_dense-p(2)).^2./(2*p(3)^2));
            y_fit = gaussdist(temp_fitvariables);
            x_cent_ideal = ideal_C2(i);
            x_cent_meas = temp_fitvariables(2);
            amp = temp_fitvariables(1);

            figure;
            plot(x,y,'b','Linewidth',2);
            hold on;
            plot(x_dense,y_fit,'k','Linewidth',2);
            plot([x_cent_ideal,x_cent_ideal],[0,amp],'--b','Linewidth',2);
            plot([x_cent_meas,x_cent_meas],[0,amp],'--k','Linewidth',2);
            ylabel('Intensity');
            xlabel('Streamwise Pixels');
            legend(["Synthetic $g_2$","Fit to $g_2$","$\mu_{2}$","$\mu_{2,m}$"],'interpreter','latex')
            grid on;    
            set(gca,'FontSize', 18);
            set(gca,'fontname','times')  % Set it to times
            xlim([60,100]);
            ylim([0,800])
        
            %Making an image for the MST journal paper
            
            pixel_offset = 26.6;

            x = x+pixel_offset;
            x_dense = x_dense+pixel_offset;
            x_cent_ideal = x_cent_ideal+pixel_offset;
            x_cent_meas = x_cent_meas+pixel_offset;

            figure;
            plot(x,y,'b-o','Linewidth',2);
            hold on;
            plot(x_dense,y_fit,'k','Linewidth',2);
            plot([x_cent_ideal,x_cent_ideal],[0,amp],':b','Linewidth',2);
            plot([x_cent_meas,x_cent_meas],[0,amp],'--k','Linewidth',2);
            ylabel('Intensity');
            xlabel('Streamwise Pixels');
            legend(["Synthetic $g_2$","Fit to $g_2$","$\mu_{i,2}$","$\mu_{i,2,m}$"],'interpreter','latex')
            grid on;    
            set(gca,'FontSize', 18);
            set(gca,'fontname','times')  % Set it to times
            xlim([102,113]);
            ylim([500,850])
        
            %with units of micrometers

            x = (x-50).*scale(1).*10^3;
            x_dense = (x_dense-50).*scale(1).*10^3;
            x_cent_ideal  = (x_cent_ideal-50).*scale(1).*10^3;
            x_cent_meas  = (x_cent_meas-50).*scale(1).*10^3;

            figure;
            plot(x,y,'b-o','Linewidth',2);
            hold on;
            plot(x_dense,y_fit,'k','Linewidth',2);
            plot([x_cent_ideal,x_cent_ideal],[0,amp],':b','Linewidth',2);
            plot([x_cent_meas,x_cent_meas],[0,amp],'--k','Linewidth',2);
            ylabel('Pixel Intensity');
            xlabel('Streamwise Distance [mm]');
            legend(["Synthetic $g_2$","Fit to $g_2$","$\mu_{i,2}$","$\mu_{i,2,m}$"],'interpreter','latex')
            grid on;    
            set(gca,'FontSize', 18);
            set(gca,'fontname','times')  % Set it to times
            xlim([2.05,2.5]);
            ylim([500,950])
        

            disp('wait');

        end


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
decay_binary = tau_fit<450;
tau_line = linspace(min(tau_fit(decay_binary)),max(tau_fit(decay_binary)),100);
tau_fit_eq = polyfit(tau_fit(decay_binary),ideal_velocity_error(decay_binary),1);
velo_line = polyval(tau_fit_eq,tau_line);
decay_error_eq = tau_fit_eq; %send this all the way to decay correction at the very end

figure;
scatter(tau_fit(decay_binary),ideal_velocity_error(decay_binary),'k','filled');
hold on;
plot(tau_line,velo_line,'r','Linewidth',2);
% title('Ideal Velocity Error vs Tau fit');
ylabel('$\Delta V_{i,\tau-err}(\tau)$ [m/s]','Interpreter','Latex')
xlabel('$\tau$ [ns]','interpreter','latex');
set(gca,'FontSize', 20);
set(gca,'fontname','times')  % Set it to times
legend(["Unique rows of data","Linear Fit"])

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