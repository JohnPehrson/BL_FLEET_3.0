function [time_av] = TimeAv_Centroid_Loading(preproc_filepath,centroids_filepath,pixel_um_resolution,runconditions_filepath,run,loc_data)
%This function will load the centroids fit from the time-averaged data from
%preprocesing. The uncorrected mean velocity is output

load(runconditions_filepath);

%% Finding out which filepath to load in
a=dir(strcat(preproc_filepath,"/*.mat"));
filenames = strings(length(a),1);
for i = 1:length(a)
filenames(i) = a(i).name;
end

%find maximum number of images processed for real data
real_filenames = filenames(contains(filenames,centroids_filepath));
real_imagenumbers = zeros(length(real_filenames),1);
for i = 1:length(real_filenames)
    localstring = convertStringsToChars(real_filenames(i));
    test = localstring((length(convertStringsToChars(centroids_filepath))+1):end-10);
    real_imagenumbers(i) = str2double(test);
end
[maxrealimages,ind] = max(real_imagenumbers);
real_filepath = strcat(preproc_filepath,'\',real_filenames(ind));

%% Loading filepath
load(real_filepath);

%% Centroids to Velocity
g1_cents = fitvariables(:,2);
g2_cents = fitvariables(:,5);
cent_times = [Delays(run,1)+Gates(run,1)/2,Delays(run,1)+Gates(run,1)+Delays(run,2)+Gates(run,2)/2];
time_dif = cent_times(2)-cent_times(1);
res = pixel_um_resolution(1);
velos = (g2_cents-g1_cents).*(res.*10^(-6))./(time_dif.*10^(-9));

%% Pixel Location to mm above the surface
time_av_velo = velos(1:loc_data(2));
time_av_heights = (length(time_av_velo)-1):-1:0;
time_av_heights = time_av_heights'.*res.*10^(-3); %mm

time_av = struct;
time_av.velo = time_av_velo;
time_av.height = time_av_heights;

%% Amplitudes
simp_amps = amplitudes(1:loc_data(2),:);
rows = size(time_av_heights,1);
tau_fit = zeros(rows,1);

%% Get gate intensities as f(y)
gate_intensity = amplitudes./Gates(run,:);

%% Decay Fitting
for i = 1:rows
    fit_ovj = fit(cent_times',gate_intensity(i,:)','exp1');
    coeffvals= coeffvalues(fit_ovj);
    tau_fit(i) = -1*inv(coeffvals(2));
end
time_av.decay = tau_fit;

end