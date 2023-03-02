clear all;clc;close all;

%% Scitech Picture maker

%% Inputs
folder_path = "C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing";
av_filepath     = "Scitech_Pictures/Time_Average_Emissions.mat";
ss_filepath     = "Scitech_Pictures/Single_Shot_Emissions.mat";
synth_filepath  = "Scitech_Pictures/Synth_Emissions.mat";


%% Data Loading
load(fullfile(folder_path,av_filepath));
av = struct;
av.x_mm = x_mm;
av.y_mm = y_mm;
av.im   = imageData_mean;

load(fullfile(folder_path,ss_filepath));
ss = struct;
ss.x_mm = x;
ss.y_mm = y;
ss.im   = imageData_ROI;

load(fullfile(folder_path,synth_filepath));
synth = struct;
synth.im   = noisey_synth_image;


%% Plotting
max_int = 3000;
fig_size = [800 100 360 800];
fontsize = 20;

f1 = figure;
f1.Position = fig_size;
    image(av.x_mm,av.y_mm,av.im);
    set(gca,'YDir','normal');
    xlabel('Streamwise Distance [mm]')
    ylabel('Height above Surface [mm]')
    axis equal;
    xlim([min(av.x_mm),max(av.x_mm)])
    ylim([min(av.y_mm),max(av.y_mm)])
    colormap(turbo(max_int));
    set(gca,'FontSize', fontsize);
    set(gca,'fontname','times')  % Set it to times
        figname = "Scitech_Pictures_Complete\TimeAverage.tif";
        saveas(gcf,figname,'tiffn')

f2 = figure;
f2.Position = fig_size;
    image(ss.x_mm,ss.y_mm,ss.im);
    set(gca,'YDir','normal');
    xlabel('Streamwise Distance [mm]')
    ylabel('Height above Surface [mm]')
    yline(0,'--w','Linewidth',2)
    axis equal;
    xlim([min(ss.x_mm),max(ss.x_mm)])
    ylim([min(ss.y_mm),max(ss.y_mm)])
    colormap(turbo(max_int));
    set(gca,'FontSize', fontsize);
    set(gca,'fontname','times')  % Set it to times
        figname = "Scitech_Pictures_Complete\SingleShot.tif";
        saveas(gcf,figname,'tiffn')

f3 = figure;
f3.Position = fig_size;
    image(ss.x_mm,ss.y_mm,synth.im);
    set(gca,'YDir','normal');
    xlabel('Synthetic Distance [mm]')
    ylabel('Synthetic Height [mm]')
    yline(0,'--w','Linewidth',2)
    axis equal;
    xlim([min(ss.x_mm),max(ss.x_mm)])
    ylim([min(ss.y_mm),max(ss.y_mm)])
    colormap(turbo(max_int));
    set(gca,'FontSize', fontsize);
    set(gca,'fontname','times')  % Set it to times
        figname = "Scitech_Pictures_Complete\Synthetic.tif";
        saveas(gcf,figname,'tiffn')




