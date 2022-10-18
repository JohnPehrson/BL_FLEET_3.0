function [dust_filter_bounds_bottom_out,dust_filter_bounds_top_out,filtered_average_ImageData_out,filtered_average_PrerunImageData] = PixelOutlierFilter(imageData_run,ROI_imagedata_flare,...
    numimages_dust_run,numimages_dust_flare, numimages_analyze,plot_switch,xlength,ylength)
%This function will identify the outlier thresholds for a run and then provide bounds
%that serve to identify outliers in individual images

    %% Loop through the prerun data, followed by the run data 
    for flare_it = 1:2 

        if flare_it ==1%flare fit/ background
            imageData = ROI_imagedata_flare;
            numimages_dust = numimages_dust_flare;
        else %run
            imageData = imageData_run;
            numimages_dust = numimages_dust_run;
        end

    
            %% Identify outliers iteratively
            numeliminated = zeros(xlength*ylength,1);
            A = 7.0;
            
            %identify rows with outliers
            pixelmean = mean(imageData,2);
            kint = (imageData-pixelmean).^4;
            k_tilda = transpose((sum(kint')/(numimages_dust-1)).^0.25);
            sigint = (imageData-pixelmean).^2;
            sigma_tilda = transpose(sqrt(sum(sigint')/(numimages_dust-1)));
            binaryoutlier = k_tilda > (sigma_tilda.*A);
            indexes = find(binaryoutlier);
            outliercount = sum(binaryoutlier);
            
            %initialize the bounds using the data set minimum/maximums
            dust_filter_bounds(:,1) = pixelmean-A.*sigma_tilda-.001;
            dust_filter_bounds(:,2) = pixelmean+A.*sigma_tilda+.001;
            
            while outliercount>0 %while there are still outliers
    
                %filter out the maximum value in the indicated rows
                for i = 1:length(indexes)
                    ind = indexes(i);
                    [maxrow,maxind] = max(imageData(ind,:));
                    
                    if plot_switch >0.5
                        pixdata = imageData(ind,:);
                        histogram(pixdata,30);
                        hold on;
                        plot([dust_filter_bounds(ind,2),dust_filter_bounds(ind,2)],[0,40]);
                        hold on;
                        plot([dust_filter_bounds(ind,1),dust_filter_bounds(ind,1)],[0,40]);
                        hold off
                        title(['Pixel Intensity Histogram - ', num2str(numimages_dust),' Pixels over a range of ',num2str(numimages_analyze),' total images']);
                        ylabel('Bin Frequency');
                        xlabel('Intensity');
                        legend('Pixels','Upper Bound','Lower Bound');
                    end
                    
                    numeliminated(indexes(i)) = numeliminated(indexes(i))+sum(imageData(indexes(i),:)>dust_filter_bounds(indexes(i),2));
                    imageData(ind,maxind) = NaN;
                    bound = (maxrow>pixelmean(ind))+1; %1 if max is less than mean (bottom bound), 2 if max is more than mean (top bound)
                    dust_filter_bounds(ind,bound) = maxrow;
                    
                end
    
                %identify rows with outliers, now considering possible NaN values
                for i = 1:size(imageData,1)
                pixelmean(i) = mean(imageData(i,~isnan(imageData(i,:))));          
                end
                kint = (imageData-pixelmean).^4;
                k_tilda = transpose((sum(kint')/(numimages_dust-1)).^0.25);
                sigint = (imageData-pixelmean).^2;
                sigma_tilda = transpose(sqrt(sum(sigint')/(numimages_dust-1)));
                binaryoutlier = k_tilda > (sigma_tilda.*A);
                indexes = find(binaryoutlier);
                outliercount = sum(binaryoutlier);
            end
            
    %         %% Plot the spatial/frequency of data eliminations
    %         numeliminated = reshape(numeliminated(:,1),[xlength,ylength]);
    %         figure;
    %         imagesc(numeliminated./numimages_dust);
    %         colormap(bone);
    %         colorbar;
    %         title('Frequency of Pixel Filtering')
    %         
            %% Dust Bounding
            dust_filter_bounds_bottom = reshape(dust_filter_bounds(:,1),[xlength,ylength]);
            dust_filter_bounds_top = reshape(dust_filter_bounds(:,2),[xlength,ylength]);
            
            %% Averaging the images and reshaping them back into one image (ignores outliers)
            binaryisreal = ~isnan(imageData);
            imageData_averaged = zeros(length(imageData(:,1)),1);
            for i = 1:length(imageData(:,1))
            imageData_averaged(i) = mean(imageData(i,binaryisreal(i,:)));
            end
            filtered_average_ImageData = reshape(imageData_averaged,[xlength,ylength]);
            
            %% Picture
            figure;
            image(filtered_average_ImageData)
            colorbar;
            colormap(turbo(4096));
            axis equal;

            %% Saving
            if flare_it ==1%flare fit/ background
                filtered_average_PrerunImageData = filtered_average_ImageData;
            else %run
                dust_filter_bounds_bottom_out = dust_filter_bounds_bottom;
                dust_filter_bounds_top_out = dust_filter_bounds_top;
                filtered_average_ImageData_out = filtered_average_ImageData;
            end


    end
end

