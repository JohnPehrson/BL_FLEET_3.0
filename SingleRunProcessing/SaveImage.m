function SaveImage(run,imageprocess_numbers,imloop,imageData_ROI,temp_centroids)
%Save Image and centroids
folder = strcat("C:\Users\clark\Documents\GitHub\BL_FLEET_2.0\RawDataProcessing\Fit_images\R",num2str(run),"\Tot_im",num2str(length(imageprocess_numbers)));
    if exist(folder, 'file') == 7 %data already exists
        %do nothing, it's already done
    else %need to create the data
        mkdir(folder);
    end
    file_num = imageprocess_numbers(imloop);
    save_filepath = strcat(folder,"\Im_num",num2str(file_num),".mat");
    save(save_filepath,"imageData_ROI","temp_centroids")
end