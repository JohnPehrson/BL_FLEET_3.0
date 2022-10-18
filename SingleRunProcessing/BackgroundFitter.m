function [background_totalfit] = BackgroundFitter(run,imageData_mean,rows,cols,numprelim_images,synth_switch)
%This function fits the background mean intensity in the data.
%This is performed in the preprocessing loop of the data

%% Expected file name
if synth_switch
savename = ['Matfiles_preprocessing/BackgroundFit_Synth_Run',num2str(run),'_',num2str(numprelim_images),'images.mat'];
else %real data
savename = ['Matfiles_preprocessing/BackgroundFit_Run',num2str(run),'_',num2str(numprelim_images),'images.mat'];
end
    %% Check if a background fit already exists
    if isfile(savename) %if it has been done before, don't redo it, just load it
        load(savename,'background_totalfit');
    else % fit doesn't already exist, compute a fit here

        %% Get the constant for the background fit
        %find the intensity of the lowest 10% of data
        sorted_pixint=sort (imageData_mean(:),1, 'ascend');
        background = mean(sorted_pixint(1:round(length(sorted_pixint)*0.10)));

%         %checking both dimensions for potential background
%         imageData_mean_colaverage = mean(imageData_mean,1);
%         imageData_mean_rowaverage = mean(imageData_mean,2);

%         figure;
%         plot(1:cols,imageData_mean_colaverage);
%         figure;
%         plot(1:rows,imageData_mean_rowaverage);
%         [XX,YY] = meshgrid(1:cols,rows:-1:1);
%         figure;
%         [C,h] = contourf(XX,YY,imageData_mean,25);
%         shading interp
%         set(h,'LineColor','none');
        
%         %% Fitting the gaussian background
%             x0 = [0,250,60,80,40];
%             LB = [0,0,30,100,x0(5)-10];
%             UB = [x0(1),rows,90,cols-100,x0(5)+20];
% 
%             %model for local noise
%             imageData_mean_backgroundsubtract = imageData_mean-background;
%             imageData_mean_backgroundsubtract(imageData_mean_backgroundsubtract<0) = 0;
%             [YY,XX] = meshgrid(1:cols,1:rows);
%             fit_gauss = @(p) (p(1)*exp(-(XX-p(2)).^2 ./ (2*p(3)^2)).*(exp(-(YY-p(4)).^2 ./ (2*p(5)^2))));
%             err_fit_gauss = @(v) fit_gauss(v)-imageData_mean_backgroundsubtract;
% 
%             %Options
%             options = optimset(@lsqcurvefit);
%             options.Display = 'off';
%             options.TolFun = 1e-5;
%             options.MaxFunEvals = 1e3;
%             options.MaxIter = 200;
%             options.FinDiffType = 'central';
% 
%             %fit
%             [localnoisecoeffs] = lsqnonlin(err_fit_gauss,x0,LB,UB,options);
%             fit_values = fit_gauss(localnoisecoeffs);    
%             streamwise_localnoise_location = localnoisecoeffs(4);
%             localnoise_subtracted = imageData_mean_backgroundsubtract-fit_values;
%             localnoise_subtracted(localnoise_subtracted<0) = 0;

%             %visualizing the fit
%             figure;
%             surf(imageData_mean);
%             title('original mean');
%             shading interp 
%             figure;
%             surf(localnoise_subtracted);
%             title('mean after background subtraction');
%             shading interp 
%             figure;
%             surf(fit_values);
%             title('background subtraction');
%             shading interp 


        %% Saving the filter for later use if needed
       background_totalfit = background; %+fit_values;
       save(savename,'background_totalfit');
    end

end