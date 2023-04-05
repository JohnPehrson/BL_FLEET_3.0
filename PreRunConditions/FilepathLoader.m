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

%% Data Saving
 currentdir  = pwd;
 idcs   = strfind(currentdir,'\');
 rootdir = currentdir(1:idcs(end)-1);
savefilepath = fullfile(rootdir,"SingleRunProcessing","TestConditions","FLEETFilePaths.mat");

save(savefilepath,'run_filepaths');
