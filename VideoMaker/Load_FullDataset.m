function [y_mm,rows_plot] = Load_FullDataset(folderpath,run,real_filename)

%% Get the folder with the most images
a=dir(strcat(folderpath,"\*"));
filenames = strings(length(a),1);
for i = 1:length(a)
filenames(i) = a(i).name;
end

%find maximum number of images processed for real data
real_filenames = filenames(contains(filenames,strcat(real_filename,num2str(run))));
real_imagenumbers = zeros(length(real_filenames),1);
real_stringlength = length(convertStringsToChars(real_filename));
for i = 1:length(real_filenames)
    localstring = convertStringsToChars(real_filenames(i));
    localstring = localstring((real_stringlength+1):end);
    real_imagenumbers(i) = str2double(localstring);
end
[maxrealimages,ind] = max(real_imagenumbers);
real_filepath = fullfile(folderpath,real_filenames(ind));

load(real_filepath);
y_mm = velocimetry_geometricloc(:,5);
allrows = 1:length(list_continuous);
rows_plot = allrows(list_continuous);

end