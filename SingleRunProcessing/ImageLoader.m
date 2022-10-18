function [imageData_out] = ImageLoader(folderName,imageName,synthfilepath,synth_switch,imagenumber)
        %% Image Loading either real or synthetic data
            if synth_switch  %synth data
            folderNameCat = strcat(synthfilepath,num2str(imagenumber,'%05.f'),".mat");
            load(folderNameCat);
            imageData_out = noisey_synth_image;
            else %real data
            folderNameCat = strcat(folderName,'\',imageName,num2str(imagenumber,'%06.f'),".tif");
            imageData_out = double(imread(folderNameCat));
            end
end