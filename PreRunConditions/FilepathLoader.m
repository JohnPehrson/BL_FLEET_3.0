clear all;close all;clc;

%% Filepath Loader
refimage_filepath = "FLEET_Filepaths.xlsx";
T = readtable(refimage_filepath);
runs = T{:,1};
folderpaths = T{:,2};
filenames = T{:,3};
folderpaths = string(folderpaths);
filenames = string(filenames);
run_filepaths = [folderpaths,filenames];

save('C:\Users\clark\Documents\GitHub\BL_FLEET_3.0\SingleRunProcessing\TestConditions/FLEETFilePaths.mat',...
    'run_filepaths');

