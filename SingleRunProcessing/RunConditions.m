clc;clear all;close all;
%% Run Conditions Script

%This script writes a mat file that contains all the run conditions for the
%3 tests. Rows denote runs

Gates = [20,700;20,700;20,700];
Delays = [100,1400;100,1400;100,1400];
Gain = [0.9;0.9;0.9];
Reynolds_approx = [6,6,6]; %million per meter

save('BLFLEETRunConditions.mat','Gates','Delays','Gain','Reynolds_approx');