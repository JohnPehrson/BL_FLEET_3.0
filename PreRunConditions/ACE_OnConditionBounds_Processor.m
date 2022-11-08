clear all;close all;clc;

%% Processing script to identify the start and end bounds for FLEET image 
%% processing based off of the ACE DAQ Reynolds number

%% Initialization of variables
 %get runs to use
    load('C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/RefData.mat');

    DAQ_times = cell(length(uniqueruns),1);
    DAQ_Res = cell(length(uniqueruns),1);
    DAQ_Machs = cell(length(uniqueruns),1);
    DAQ_Velos = cell(length(uniqueruns),1);
    DAQ_start_stops = zeros(length(uniqueruns),2);
    DAQ_start_stops(:,2) = [24.4;26.4;24.3;28.1;26.9;...
                            27.6;26.6;28.5;19.5;30.4;25.9;...
                            25.9;24.8;26.7];
    DAQ_flare_start_stop = zeros(length(uniqueruns),2);
    target_Re = 6.1*10^6;
    target_Re_unc = 0.1*10^6;
    Run_Re = zeros(length(uniqueruns),2); %left column is mean, right column is 2 std
    Run_Mach = zeros(length(uniqueruns),2); %left column is mean, right column is 2 std
    Run_Mean_Velo = zeros(length(uniqueruns),2); %left column is mean, right column is 2 std


%only load ACE DAQ in if I haven't before
savename_DAQ = "DAQ_Reynolds.mat";
    if isfile(savename_DAQ) %I already have the data, just load the file in
        load(savename_DAQ);
    else % Need to load in the data
    
    %% Load in Reynolds numbers and time from each run
        %get the filepaths for the ACE DAQ files with those names
        ACE_generalfilepath = 'D:/ACEDAQ';
            %get all filenames of runs
            a=dir(strcat(ACE_generalfilepath, "/*.xlsx"));
            all_filenames = strings(length(uniqueruns),1);
            for i = 1:length(a)
                all_filenames(i) = a(i).name;
            end
            %get only the appropriate runs
            used_filepaths = strings(length(uniqueruns),1);
            for j = 1:length(uniqueruns)
                which_run_binary = contains(all_filenames,num2str(uniqueruns(j)));
                used_filepaths(j) = all_filenames(which_run_binary);
            end
            %Load in those runs and save out the appropriate data
            for k = 1:length(used_filepaths)
                full_filepath = strcat(ACE_generalfilepath,'/',used_filepaths(k));
                single_DAQ_time = readmatrix(full_filepath,'Sheet','Reduced Data','Range','A2:A1000'); %seconds
                single_DAQ_Re = readmatrix(full_filepath,'Sheet','Reduced Data','Range','AY2:AY1000');
                single_DAQ_Mach = readmatrix(full_filepath,'Sheet','Reduced Data','Range','AZ2:AZ1000');
                single_DAQ_Velos = readmatrix(full_filepath,'Sheet','Reduced Data','Range','AW2:AW1000');
    
                %put the single run data into a macro-variable
                DAQ_times{k} = single_DAQ_time;
                DAQ_Res{k} = single_DAQ_Re;
                DAQ_Machs{k} = single_DAQ_Mach;
                DAQ_Velos{k} = single_DAQ_Velos;
            end

            save(savename_DAQ,'DAQ_times','DAQ_Res','DAQ_Machs','DAQ_Velos');
    end %end of the if loop to load in data

%% Finding the steady state start and stop time for the ACE Tunnel
    %focus on visualization to make sure it appears reasonable

    %Find start time as the first time the Reynolds passes 6M/m
    for i = 1:length(uniqueruns)
        single_DAQ_Re = DAQ_Res{i};
        single_DAQ_time = DAQ_times{i};
        above_condition = single_DAQ_Re>(target_Re-target_Re_unc);
        V1 = find(above_condition, 1, 'first');
        DAQ_start_stops(i,1) = single_DAQ_time(V1);
    end
    
    %Find the end time based on a sufficiently strong Reynolds gradient
    Re_grad = cell(length(uniqueruns),1);
    LocalMaxes = cell(length(uniqueruns),1);
    for i = 1:length(uniqueruns)
        single_DAQ_Re = DAQ_Res{i};
        single_DAQ_time = DAQ_times{i};
            single_Re_grad = gradient(single_DAQ_Re);
            Re_grad{i} = single_Re_grad;
            TF_max = islocalmax(single_Re_grad,'MinProminence',2.5*10^5);
            TF_min = islocalmin(single_Re_grad,'MinProminence',2*10^5);
            LocalMaxes{i} = or(TF_max,TF_min);
    end

    %Find the end time for the pre-run section via the reynolds gradient
    for i = 1:length(uniqueruns)
        single_DAQ_Re = DAQ_Res{i};
        single_DAQ_time = DAQ_times{i};
        startup = find(LocalMaxes{i}, 1, 'first');
        pre_startup = startup-3;
        DAQ_flare_start_stop(i,2) = single_DAQ_time(pre_startup)*1000;
    end

    %Plotting
    figure(1);
    for i = 1:length(uniqueruns)
        subplot(2,1,1);
        single_times = DAQ_times{i};
        single_Res = DAQ_Res{i};
            plot(single_times,single_Res,'k','Linewidth',2);
            hold on;
            %condition visualization
            plot(DAQ_times{i},target_Re.*ones(length(DAQ_times{i}),1),'b','Linewidth',1)
            plot(DAQ_times{i},(target_Re-target_Re_unc).*ones(length(DAQ_times{i}),1),'--b','Linewidth',1)
            plot(DAQ_times{i},(target_Re+target_Re_unc).*ones(length(DAQ_times{i}),1),'--b','Linewidth',1)

            %start time
            plot([DAQ_start_stops(i,1),DAQ_start_stops(i,1)],[0,8*10^6],'g','Linewidth',2)
            %end time
            plot([DAQ_start_stops(i,2),DAQ_start_stops(i,2)],[0,8*10^6],'r','Linewidth',2)         
            hold off;

        subplot(2,1,2);
        single_Re_grad = Re_grad{i};
        plot(single_times,single_Re_grad,'k','Linewidth',2);
        hold on;
        scatter(single_times(LocalMaxes{i}),single_Re_grad(LocalMaxes{i}),'r*')
        hold off;

    end
close all;

%% Report the Reynolds number and Mach numberwith uncertainty bounds for each run

for i = 1:length(uniqueruns)
    single_times = DAQ_times{i};
    single_Res = DAQ_Res{i};
    single_Machs = DAQ_Machs{i};
    single_Velos = DAQ_Velos{i};

    binary_run = (single_times>=DAQ_start_stops(i,1))&(single_times<=DAQ_start_stops(i,2));
    mean_Re = mean(single_Res(binary_run));
    mean_Mach = mean(single_Machs(binary_run));
    mean_Velos = mean(single_Velos(binary_run));

    doublestd_Re = 2*std(single_Res(binary_run));
    doublestd_Mach = 2*std(single_Machs(binary_run));
    doublestd_Velos = 2*std(single_Velos(binary_run));

    Run_Re(i,:) = [mean_Re,doublestd_Re];
    Run_Mach(i,:) = [mean_Mach,doublestd_Mach];
    Run_Mean_Velo(i,:) = [mean_Velos,doublestd_Velos];
end

Run_Re_scaled = Run_Re/(10^6);

%% Plotting individual runs 'steady state' Mach and Reynolds
for i = 1:length(uniqueruns)
    single_times = DAQ_times{i};
    single_Res = DAQ_Res{i};
    single_Machs = DAQ_Machs{i};

    binary_run = (single_times>=DAQ_start_stops(i,1))&(single_times<=DAQ_start_stops(i,2));
    single_times = single_times(binary_run);
    single_Res = single_Res(binary_run);
    single_Res = single_Res/(10^6);
    single_Machs = single_Machs(binary_run);

    figure(2);
    colororder({'k','r'})
    yyaxis left
    plot(single_times,single_Machs,'k','Linewidth',2);
    ylim([5.5,6]);
    ylabel('Mach Number [-]')
    xlabel('Time since tunnel start [s]')
    title(['Run',num2str(uniqueruns(i))]);
    grid on;
    set(gca,'FontSize', 17);
    set(gca,'fontname','times')  % Set it to times
    
    yyaxis right
    plot(single_times,single_Res,'r','Linewidth',2);
    ylim([5.5,7.5]);
    ylabel('Reynolds Number [10^6/m]')
    grid on;
    set(gca,'FontSize', 17);
    set(gca,'fontname','times')  % Set it to times

    legend('Mach','Reynolds')
end

%scale start and end time to 1000fps
DAQ_start_stops = DAQ_start_stops.*1000;

%start with the first image for the flare
DAQ_flare_start_stop(:,1) = 1;



%% Save out Data
save('C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/ACE_Data.mat',...
    'Run_Re_scaled','DAQ_start_stops','DAQ_flare_start_stop','Run_Mach','Run_Mean_Velo');

