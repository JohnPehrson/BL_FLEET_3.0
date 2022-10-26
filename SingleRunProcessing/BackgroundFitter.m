function [background_totalfit] = BackgroundFitter(run,imageData_mean,rows,cols,numprelim_images,synth_switch)
%This function fits the background mean intensity in the data.
%This is performed in the preprocessing loop of the data

%% Expected file name
if synth_switch
savename = ['Preprocessing_Filestorage/BackgroundFit_Synth_Run',num2str(run),'_',num2str(numprelim_images),'images.mat'];
else %real data
savename = ['Preprocessing_Filestorage/BackgroundFit_Run',num2str(run),'_',num2str(numprelim_images),'images.mat'];
end
    %% Check if a background fit already exists
    if isfile(savename) %if it has been done before, don't redo it, just load it
        load(savename,'background_totalfit');
    else % fit doesn't already exist, compute a fit here

        %% Get the constant for the background fit
        %find the intensity of the lowest 10% of data
        sorted_pixint=sort (imageData_mean(:),1, 'ascend');
        background = mean(sorted_pixint(1:round(length(sorted_pixint)*0.10)));

        %% Saving the filter for later use if needed
       background_totalfit = background; %+fit_values;
       save(savename,'background_totalfit');
    end

end