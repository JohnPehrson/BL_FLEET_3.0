function [image_data_all,centroids_all,num_ims] = LoadImages(real_filepath)
    %Loading images
    ims_folderpath = strcat(real_filepath);
    a=dir(strcat(ims_folderpath,"\*"));
    filenames = strings(length(a),1);
    for i = 3:(length(a))
    filenames(i) = a(i).name;
    end
    filenames = filenames(3:end);
    
    num_ims = length(filenames);
    
    for i = 1:length(filenames)
        if i==1
            load(fullfile(real_filepath,filenames(i)));
            image_data_all = zeros([size(imageData_ROI),length(filenames)]);
            centroids_all = zeros([size(temp_centroids),length(filenames)]);
            image_data_all(:,:,i) = imageData_ROI;
            centroids_all(:,:,i) = temp_centroids;
        else
            load(fullfile(real_filepath,filenames(i)));
            image_data_all(:,:,i) = imageData_ROI;
            centroids_all(:,:,i) = temp_centroids;
        end
    
    end
end