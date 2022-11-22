function [fullfilepath,run_start_end,nondim_velo_error,synth_input_tau_fit,...
    synth_input_velocity_mean] = Synthetic_Data_Gen(runconditions_filepath,...
    resolution_filepath,run,synth_real_replicate,synth_numimages,ROI,fitting_limits)
%This function will generate a set of synthetic data if one does not
%already exist that matches the desired parameters. 

%% Load in real data 
load(runconditions_filepath);
load(resolution_filepath);
load(synth_real_replicate);
rows = transpose(1:length(velocity_mean));


%% Set up folder for synthetic data saving
folder = 'C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\SynthData';
subfolder = ['Run',num2str(run),'_SynthImages',num2str(synth_numimages)];

%check to see if the folder already exists, but is empty
if exist(fullfile(folder, subfolder), 'file') == 7 %data already exists
        a=dir(fullfile(folder, subfolder));
        out=size(a,1);
    if out==(synth_numimages+2)%enough files in folder
        %don't reprocess
    else %delete and reprocess
        %delete files in folder
        for i= 3:length(a)
            filename = a(i).name;
            delete(fullfile(folder, subfolder,filename));
        end

        %delete folder
        rmdir(fullfile(folder, subfolder));    
    end
end

    if exist(fullfile(folder, subfolder), 'file') == 7 %data already exists
        %do nothing, it's already done
    else %need to create the data
        mkdir(folder,subfolder);

                %Generate the ideal synthetic data once
                [ideal_synth_image,nondim_velo_error,synth_input_tau_fit,synth_input_velocity_mean,...
                ROI_mean_SNR,ROI_gate1_location_bounds,ROI_gate2_location_bounds] = SynthData_Ideal(rows,...
                                        velocity_mean,velocimetry_geometricloc(:,5),pixel_um_resolution(run,:),...
                                        Gates(run,:),Delays(run,:),doublegauss_fitvariables,...
                                        emissionlocatingdata,tau_fit,ROI,gate1_location_bounds,...
                                        gate2_location_bounds,fitting_limits,list_continuous,mean_SNR,...
                                        run);

            for i = 1:synth_numimages
                %add in noise for every image
                [noisey_synth_image] = SynthData_NoiseApplication(ideal_synth_image,ROI_mean_SNR,ROI_gate1_location_bounds,...
                                        ROI_gate2_location_bounds,background_totalfit,doublegauss_fitvariables(:,4),list_continuous);
                
                filenumber = sprintf( '%05d', i ) ;
                filename = ['Synth_Data_',num2str(filenumber)];
                fullfilepath = [folder,'\',subfolder,'\',filename];
                save(fullfilepath,'noisey_synth_image','nondim_velo_error','synth_input_tau_fit','synth_input_velocity_mean')
            end
    end
    fullfilepath = [folder,'\',subfolder,'\','Synth_Data_'];
   
    if exist(fullfile(folder, subfolder), 'file') == 7 %data already exists
        singlefilepath = [fullfilepath,'00001.mat'];
    load(singlefilepath,'noisey_synth_image','nondim_velo_error','synth_input_tau_fit','synth_input_velocity_mean');
    end

    run_start_end(run,:) = [1,synth_numimages];

end

