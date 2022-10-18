function [real,synth,maxrealimages] = RealSynthDataLoader(folderpath,real_filename,real_stringlength,...
                                    synth_filename,synth_stringlength)
%This function identifies the real and synthetic data with the most images
%processed and then loads them into memory

%% find file paths
a=dir(strcat(folderpath, "/*.mat"));
filenames = strings(length(a),1);
for i = 1:length(a)
filenames(i) = a(i).name;
end

%find maximum number of images processed for real data
real_filenames = filenames(contains(filenames,real_filename));
real_imagenumbers = zeros(length(real_filenames),1);
for i = 1:length(real_filenames)
    localstring = convertStringsToChars(real_filenames(i));
    localstring = localstring((real_stringlength+1):end-4);
    real_imagenumbers(i) = str2double(localstring);
end
[maxrealimages,ind] = max(real_imagenumbers);
real_filepath = strcat(folderpath,'\',real_filenames(ind));

%find maximum number of images processed for synthetic data
synth_filenames = filenames(contains(filenames,synth_filename));
synth_imagenumbers = zeros(length(synth_filenames),1);
for i = 1:length(synth_filenames)
    localstring = convertStringsToChars(synth_filenames(i));
    localstring = localstring((synth_stringlength+1):end-4);
    synth_imagenumbers(i) = str2double(localstring);
end
[im,ind] = max(synth_imagenumbers);
synth_filepath = strcat(folderpath,'\',synth_filenames(ind));

%% Load in data
%real
real = load(real_filepath);
%synthetic
synth = load(synth_filepath);
end