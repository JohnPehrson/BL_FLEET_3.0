clc;clear all;close all;

%This script writes a mat file that contains all the run conditions for the
%tests, different rows correspond to different tests

%% Run Conditions
Gates = [25,700;...
         20,700;...
         20,700;...
         20,700;...
         15,800;...
         12,800;...
         10,800;...
         10,800;...
         10,800;...
         10,800;...
         10,800;...
         10,800;...
         10,800;...
         10,800];

Delays = [40,1400;...
          40,1400;...
          40,1400;
          40,1400;
          40,1600;
          40,1600;
          50,1600;
          50,1500;
          50,1500;
          50,1500;
          50,1500;
          50,1500;
          50,1500;
          50,1500];

%positive is upstream, negative is downstream, second column is uncertainty
SpanwiseCameraAngle = [1.5,1;
                       1.5,1;
                       2.5,1;
                       2.5,1;
                       1,1;
                       0,1;
                       0,1;
                       -3,1;
                       -2.5,1;
                       0,1;
                       0,1;
                       0,1;
                       -2.5,1;
                       0,1;]; 

RunNames =  ["SRA_22C";...
            "SRA_22C";...
            "PB_22C";...
            "PB_22C";...
            "PB_22C";...
            "PB_32C";...
            "PB_32R";...
            "PB_33C";...
            "PB_21C";...
            "PB_11C";...
            "SRA_11C";...
            "SRA_11C";...
            "SRA_21C";...
            "SRA_22C"];

%% Location Stuff
hole_x = [181;220;220;320;320;320;...
          370;370;370;410;410;410;...
          510;510;510;540;...
          540;570;570;570];
trip_d = 6.84;
hole_y = [0;0;17.1;0;17.1;34.2;...
          0;-46;46;-38;0;38;...
          0;-7.75;7.75;0;...
          7.75;0;7.75;-7.75];

%hole labels
holelabels = ["11C";"12C";"12R";"21C";"21R";"21RR";"22C";"22L";"22R";...
               "23L";"23C";"23R";"31C";"31L";"31R";"32C";"32R";...
               "33C";"33R";"33L";];

downstream_loc = zeros(length(RunNames),1);
spanwise_loc = zeros(length(RunNames),1);
downstream_loc_unc = 1.*ones(length(RunNames),1);
spanwise_loc_unc = 1.*ones(length(RunNames),1);

for i = 1:length(holelabels)
    singlehole = holelabels(i);
    TF = contains(RunNames,singlehole);
    downstream_loc(TF) = hole_x(i);
    spanwise_loc(TF) = hole_y(i);
end

%% Data Saving
 currentdir  = pwd;
 idcs   = strfind(currentdir,'\');
 rootdir = currentdir(1:idcs(end)-1);
savefilepath = fullfile(rootdir,"SingleRunProcessing","TestConditions","BLFLEETRunConditions.mat");

save(savefilepath,'Gates','Delays','RunNames','downstream_loc','spanwise_loc','downstream_loc_unc','spanwise_loc_unc',...
    'SpanwiseCameraAngle');




