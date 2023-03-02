function [aw,no_aw] = LoadCFD_Data_Sept2022(inlet_cfd_filepath)
%This function loads in and sorts the RANS CFD. CFD was run with the
%SARCQCR RANS model with both adiabatic wall and non-adiabtic wall
%assumption
    
    no_aw = struct;
    aw = struct;
    no_aw.u = zeros(400,6);
    no_aw.h = zeros(400,6);
    aw.u = zeros(400,6);
    aw.h = zeros(400,6);
    
    no_aw_all = readmatrix(inlet_cfd_filepath,'Sheet','SARCQCR','Range','A4:Q404'); %seconds
    no_aw.u = no_aw_all(:,1:3:end);
    no_aw.h = no_aw_all(:,2:3:end);
    
    aw_all = readmatrix(inlet_cfd_filepath,'Sheet','SARCQCR','Range','S4:AI404'); %seconds
    aw.u = aw_all(:,1:3:end);
    aw.h = aw_all(:,2:3:end);

end