function [tau_fit] = IntensityDecay(imageData_mean,time_averaged_fit,scale,...
    gate1_location_bounds,gate2_location_bounds,Gates,Delays,fitting_limits,nearwall_bounds,...
    background_totalfit,amplitudes)
%This identifies the decay constant as a function of height using the
%time averaged FLEET emissions

%% Initialize variables
t = [ Delays(1)+Gates(1)/2 ,Delays(1)+Gates(1)+Delays(2)+Gates(2)/2]; %ns
rows = size(imageData_mean,1);
cols = size(imageData_mean,2);
gate_intensity = zeros(rows,2);
fitdata = imageData_mean-time_averaged_fit;
notwall_bounds = [1:(nearwall_bounds(1)-1)];
closewall_bounds = nearwall_bounds(1):nearwall_bounds(2);
fitvariables = zeros(rows,6);
tau_fit = zeros(rows,1);

%% Get gate intensities as f(y)
    gate_intensity = amplitudes./Gates;

%% Decay Fitting
for i = 1:rows
    fit_ovj = fit(t',gate_intensity(i,:)','exp1');
    coeffvals= coeffvalues(fit_ovj);
    tau_fit(i) = -1*inv(coeffvals(2));
end

%     %% Visualizing the fit
%     yplot = 1:rows;
%     figure;
%     plot(tau_fit(yplot),yplot);
%     hold on;
%     xlabel('Intensity Decay Constant');
%     ylabel('Pixels in the Image');
%     title('Decay Constant Fitting Real Data');
%     set(gca, 'YDir','reverse')
%     

end

