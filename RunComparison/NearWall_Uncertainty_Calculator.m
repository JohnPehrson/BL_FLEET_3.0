function [fitobject] = NearWall_Uncertainty_Calculator(near_wall_folder_path,near_wall_file_name,Gates,Delays,pixel_um_resolution)
%This function loads in data, calculates velocity, and then fits the
%various velocity profiles (of the same run with different near-wall fit
%parameters) with a 95% prediction interval to show the impact of near-wall
%effects

    %% Data Loading 
        a=dir(strcat(near_wall_folder_path,near_wall_file_name, "*.mat"));
        all_filenames = strings(size(a,1),1);
        for i = 1:length(a)
            all_filenames(i) = a(i).name;
        end
        near_wall_filename = strings(size(a,1),1);
        for i = 1:length(near_wall_filename)
            rr_name = strcat(near_wall_file_name,num2str(i),"_");
            bin_test = contains(all_filenames,rr_name);
            near_wall_filename(i) = all_filenames(bin_test);
        end

     %% near-wall stuff
            nw_run = convertStringsToChars(near_wall_filename(1));
            nw_run2 = str2num(convertCharsToStrings(nw_run(25)));
            nw_centroids        = cell(length(near_wall_filename),1);
            nw_centroid_erros   = cell(length(near_wall_filename),1);
            nw_heights          = cell(length(near_wall_filename),1);
            nw_gates            = cell(length(near_wall_filename),1);
            nw_delays           = cell(length(near_wall_filename),1);
            nw_resolutions      = cell(length(near_wall_filename),1);

            nw_velocities                   = cell(length(near_wall_filename),1);
            nw_random_velocity_error        = cell(length(near_wall_filename),1);
            nw_systematic_velocity_error    = cell(length(near_wall_filename),1);

            nw_allvelocities    = [];
            nw_allheights       = [];

    %% Loading in Data from the Sept 2022 Campaign, same run with different fit parameters near the wall
        for i= 1:length(near_wall_filename)
            load(fullfile(near_wall_folder_path,near_wall_filename(i)))
            nw_centroids{i} = centroids;
            nw_centroid_erros{i} = centroid_error;
            nw_heights{i} = y_mm;
            nw_gates{i} = Gates(nw_run2,:);
            nw_delays{i} = Delays(nw_run2,:);
            nw_resolutions{i} = pixel_um_resolution(nw_run2,:);
        end

    % Calculating Velocity
        for i = 1:length(near_wall_filename)
            run_cent        = nw_centroids{i};
            run_cent_error  = nw_centroid_erros{i};
            run_resolution  = nw_resolutions{i};
            run_gates       = nw_gates{i};
            run_delays      = nw_delays{i};
            y_mm            = nw_heights{i};
        
        
            [velocity,random_velocity_error,systematic_velocity_error] = VelocityFinder(run_cent,...
                run_cent_error,run_resolution,run_gates,run_delays);
        
        nw_velocities{i}                 = velocity;
        nw_random_velocity_error{i}      = random_velocity_error;
        nw_systematic_velocity_error{i}  = systematic_velocity_error;
        nw_allvelocities = [nw_allvelocities,velocity];
        nw_allheights = [nw_allheights,(y_mm')];%+eps*i*100)];
        end

    %% Prediction Interval 
        predict_interval = zeros(size(nw_allheights,1),3);
        for i= 1:size(nw_allheights,1)
            height = nw_allheights(i,1);
            velos = nw_allvelocities(i,:);
            mean_velo = mean(velos);

            n = length(velos);
            SSE = sum((velos-mean_velo).^2);
            MSE = SSE/(n-2);
            s2_predict = MSE*(1+1/n);
                p = .95;   
                t = tinv(p,n-2);
            side = t*sqrt(s2_predict);
            predict_interval(i,:)= [mean_velo-side,mean_velo,mean_velo+side];
        end

    %% Calculating the uncertainty as a function of height
    unc_heights = nw_allheights(:,1);
    near_wall_unc = predict_interval(:,2)-predict_interval(:,1);

    %% Make a function to report the uncertainty in the mean velocity as a function of height (so I can apply this same thing to other runs)
    unc_heights_aw = unc_heights(unc_heights>=0);
    near_wall_unc_aw = near_wall_unc(unc_heights>=0);
    fitobject = fit(unc_heights_aw,near_wall_unc_aw,'cubicinterp');
    hypo_heights = 0:0.01:12;


%     figure;
%     plot(unc_heights_aw,near_wall_unc_aw);
%     hold on;
%     plot(fitobject,unc_heights_aw,near_wall_unc_aw)
%     hold on;
%     hypo_error = feval(fitobject,hypo_heights);
%     plot(hypo_heights,hypo_error)



%     %% Plotting
%         figure;
%         plot(nw_allvelocities,nw_allheights);
%         hold on;
%         plot(predict_interval(:,2),nw_allheights(:,1),'k','Linewidth',1.5);
%         plot(predict_interval(:,1),nw_allheights(:,1),'--k','Linewidth',1.5);
%         plot(predict_interval(:,3),nw_allheights(:,1),'--k','Linewidth',1.5);
%         grid on;    
%         title('Near Wall Uncertainty Prediction Interval');
%         xlabel('Velocity [m/s]');
%         ylabel('Height above the surface [mm]');
%         set(gca,'FontSize', 15);
%         set(gca,'fontname','times')  % Set it to times
%         xlim([0,600]);
%         ylim([0,1])
%        

end