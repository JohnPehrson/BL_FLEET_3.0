function [filt_centroids] = Load_filtered_centroids(folderpath,run,filt_filename)
    folder_combined = strcat(folderpath,num2str(run));
    a=dir(strcat(folder_combined,"\*"));
    filenames = strings(length(a),1);
    for i = 1:length(a)
    filenames(i) = a(i).name;
    end
    
    %find maximum number of images processed for real data
    real_filenames = filenames(contains(filenames,filt_filename));
    real_imagenumbers = zeros(length(real_filenames),1);
    real_stringlength = length(convertStringsToChars(filt_filename));
    for i = 1:length(real_filenames)
        localstring = convertStringsToChars(real_filenames(i));
        localstring = localstring((real_stringlength+1):end-4);
        real_imagenumbers(i) = str2double(localstring);
    end
    [~,ind] = max(real_imagenumbers);
    real_filepath = strcat(folderpath,num2str(run),'\',real_filenames(ind));
    load(real_filepath); %has the output in here

end