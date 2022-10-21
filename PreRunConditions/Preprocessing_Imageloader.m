function [imageData_out] = Preprocessing_Imageloader(generalfoldername,runfolder,imageName,imagenumber)
        %% Image Loading either real or synthetic data
            folderNameCat = strcat(generalfoldername,'\',imageName,'\',runfolder,num2str(imagenumber,'%06.f'),".tif");
            imageData_out = double(imread(folderNameCat));

end