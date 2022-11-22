function [ideal_image_ROI,gate2_trimmed,gate1_trimmed,gate2_image_ROI,...
         ROI_mean_SNR,ROI_gate1_location_bounds,ROI_gate2_location_bounds,...
         loc_midtime,Amp_g2_mod] = SynthGateGen(p_w,...
         Amp_g1,Width_g1,C1,Amp_g2,Width_g2,...
         count_cols,count_rows,emissionlocatingdata,...
         Delays,Gates,tau_fit,...
         velocity_mean,scale,ROI,mean_SNR,...
         gate1_location_bounds,gate2_location_bounds)
%Function to generate the IRO gates based on amplitude decay, diffusion (width), 
%and displacement (centroid). This happens as a function of height above the surface, where values change.
%This is meant to be used to model FLEEET emissions

%Clark Pehrson
%Mar 17, 2022

%% Definitions
timesteps = (Delays(1)+Gates(1)+Delays(2)+Gates(2)); %ns
t = [Delays(1)+Gates(1)/2 ,Delays(1)+Delays(2)+Gates(1)+Gates(2)/2];
time_stepper = 1:timesteps;
Amp_g1_mod = smooth(Amp_g1,6,'lowess');
Amp_g2_mod = smooth(Amp_g2,6,'lowess');
mean_SNR = smooth(mean_SNR,6,'lowess');
Amp_gates = [Amp_g1_mod,Amp_g2_mod];


%% Increasing gate amplitudes to counteract the effect of gate spotiness
% 
% %gate 1
% [spotiness_percent_g1] = GateSpotiness(1,0);
% Amp_gates(:,1) = Amp_g1_mod./spotiness_percent_g1;
% 
% %gate 2
% [spotiness_percent_g2] = GateSpotiness(2,mean_SNR);
% Amp_gates(:,2) = Amp_g2_mod./(spotiness_percent_g2);

%% Change the decay constant in simulation based on the modifications to the amplitude
tau_gen = zeros(size(tau_fit));
for i = 1:length(tau_fit)
test_exp = fit(t',Amp_gates(i,:)','exp1');
coeffvals= coeffvalues(test_exp);
tau_gen(i) = -1.5/coeffvals(2);
end
tau_gen(end) = tau_gen(end-1);

%% Fundamental Emission Shape
emissionshape = @(amp_g,loc_g,wid_g,x) (amp_g*exp(-(x-loc_g).^2 ./ (2*wid_g^2)));

%% Time-based Equations
%amp
amp_func = @(a_0,t_1,time) a_0.*exp(-time./(t_1)); %single exponential function
width_func = @(wid,time) wid*polyval(p_w,time*10^-9);

%% Initializing Time-based Matraxies (rows are rows in the image, columns are time steps in the simulation)
%IRO and amplitude
IRO_amp_timed = zeros(count_rows,timesteps);
for i = 1:count_rows %for each row
    for j = 1:length(Gates) %for both gates in a row
            amplitude_decay = amp_func(Amp_gates(i,j),tau_gen(i),1:Gates(j));
            if j==1 %gate 1
                gate_time = (Delays(j)):(Delays(j)+Gates(j)-1);
            else  %gate 2
                gate_time = (Delays(1)+Gates(1)+Delays(j)):(Delays(1)+Gates(1)+Delays(j)+Gates(j)-1);
            end
            IRO_amp_timed(i,gate_time) = amplitude_decay;
    end
end
binary_isData = IRO_amp_timed>0;

%locations due to velocity
loc_mat_timed = zeros(count_rows,timesteps);
rel_time = (Delays(1)+Gates(1)/2);
for i = 1:count_rows %for each row
    row_vel = velocity_mean(i);
    loc_mat_timed(i,:) = (time_stepper-rel_time).*row_vel./(scale(1)*10^9);
end
loc_mat_timed = loc_mat_timed+C1;
loc_mat_timed(~binary_isData) = 0;
loc_midtime = loc_mat_timed(:,t(2)); %the time at the time-based midpoint of the gate


%width due to velocity
width_mat_timed = zeros(count_rows,timesteps);
rel_time = (Delays(1)+Gates(1)/2);
for i = 1:count_rows %for each row
    width_mat_timed(i,:) = Width_g1(i)*polyval(p_w,(time_stepper-rel_time)*10^-9);
end
width_mat_timed(~binary_isData) = 0;

%% Superposition of individual time-stepped elements
xloc = 1:count_cols;
x_pix_int = zeros(1,length(xloc));

hasdata_time = time_stepper(binary_isData(1,:));
ideal_image = zeros(count_rows,count_cols);
gate1_trimmed = zeros(count_rows,count_cols);
gate2_trimmed = zeros(count_rows,count_cols);
for j = 1:count_rows
    gauss_timeindep = zeros(length(hasdata_time),length(xloc));
    for i = 1:length(hasdata_time)
    gauss_timeindep(i,:) = emissionshape(IRO_amp_timed(j,hasdata_time(i)),loc_mat_timed(j,hasdata_time(i)),width_mat_timed(j,hasdata_time(i)),xloc);
    end
    ideal_image(j,:) = sum(gauss_timeindep,1); %superposition
    gate1_trimmed(j,:) = sum(gauss_timeindep(1:Gates(1),:),1); %superposition
    gate2_trimmed(j,:) = sum(gauss_timeindep((Gates(1)+1):end,:),1); %superposition
end

%% Peter-out of the intensity
row_peter = 50;
binary_amp = linspace(0.20,1,row_peter);
ideal_image(1:row_peter,:) = ideal_image(1:row_peter,:).*binary_amp';

% plotting
figure;
image(ideal_image)
colorbar;
colormap(bone(round(max(ideal_image(:)))));
set(gca,'FontSize', 15);
set(gca,'fontname','times')  % Set it to times
title('Synthetic Replication of Reported Data');

%% Resize the image to be in a hypothetical ROI
ideal_image_ROI             = zeros(length(ROI(1):ROI(2)),length(ROI(3):ROI(4)));
gate2_image_ROI             = zeros(length(ROI(1):ROI(2)),length(ROI(3):ROI(4)));
ROI_mean_SNR                = zeros(length(ROI(1):ROI(2)),1);
ROI_gate1_location_bounds   = zeros(length(ROI(1):ROI(2)),3);
ROI_gate2_location_bounds   = zeros(length(ROI(1):ROI(2)),3);

cols_ROI = size(ideal_image_ROI,2);
rows_ROI = size(ideal_image_ROI,1);
count_rows = size(ideal_image,1);

%place the synth data in the ROI using the known emission locating data
trimrows = count_rows;
rows_offset = ceil(emissionlocatingdata(2)-trimrows);
rows_inset = rows_offset:(rows_offset+trimrows-1);
ideal_image_ROI(rows_inset,:)           = ideal_image(1:trimrows,:);
gate2_image_ROI(rows_inset,:)           = gate2_trimmed(1:trimrows,:);
ROI_mean_SNR(rows_inset,:)              = mean_SNR(1:trimrows,:);
ROI_gate1_location_bounds(rows_inset,:) = gate1_location_bounds(1:trimrows,:);
ROI_gate2_location_bounds(rows_inset,:) = gate2_location_bounds(1:trimrows,:);

%repeat with the upside-down replica of the image
ideal_image_rev                 = flipud(ideal_image(1:trimrows,:));
ideal_g2_rev                    = flipud(gate2_trimmed(1:trimrows,:));
ROI_mean_SNR_rev                = flipud(mean_SNR(1:trimrows,:));
ROI_gate1_location_bounds_rev   = flipud(gate1_location_bounds(1:trimrows,:));
ROI_gate2_location_bounds_rev   = flipud(gate2_location_bounds(1:trimrows,:));

nextrow = (rows_inset(end));
    %take the less of the two lengths
    max_fill = min([size(ideal_image_ROI,1),nextrow+trimrows]);
fill_rows = nextrow:(max_fill-1);
ideal_image_ROI(fill_rows,:)            = ideal_image_rev(1:length(fill_rows),:);
gate2_image_ROI(fill_rows,:)            = ideal_g2_rev(1:length(fill_rows),:);
ROI_mean_SNR(fill_rows,:)               = ROI_mean_SNR_rev(1:length(fill_rows),:);
ROI_gate1_location_bounds(fill_rows,:)  = ROI_gate1_location_bounds_rev(1:length(fill_rows),:);
ROI_gate2_location_bounds(fill_rows,:)  = ROI_gate2_location_bounds_rev(1:length(fill_rows),:);

end_rows = [1:rows_inset(1),fill_rows(end):length(ROI_mean_SNR)];
for i = 1:length(end_rows)
    ROI_gate1_location_bounds(end_rows(i),:) = ROI_gate1_location_bounds(rows_inset(1),:);
    ROI_gate2_location_bounds(end_rows(i),:) = ROI_gate2_location_bounds(rows_inset(1),:);
end

figure;
image(ideal_image_ROI)
colorbar;
colormap(bone(round(max(ideal_image_ROI(:)))));
set(gca,'FontSize', 15);
set(gca,'fontname','times')  % Set it to times
title('Synthetic Replication of Image ROI');

end

