function [cfd_turb_prof] = LoadCFD_turb_BL(CFD_turbulent_profile_filepath,resolution)
%This function loads in the turbulent boundary layer from CFD

cfd_turb_prof = struct;

T = readtable(CFD_turbulent_profile_filepath);
y = T{:,1}.*1000;
v = T{:,2};

cfd_turb_prof.y_pix = y.*resolution./2;
cfd_turb_prof.v_nondim = v./max(v);

end