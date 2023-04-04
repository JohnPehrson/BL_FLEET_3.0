function [norm_velo_defect,y_norm,norm_velo_defect_max,norm_velo_defect_min] = ...
    OuterVariables(u_tau_c,MeanVelo_Profile,MeanVelo_Profile_uncertainty,EdgeVelo,height,BL_height)
%This function scales boundary layer profiles to the outer variables

%velocity defect
v_star = u_tau_c;
norm_velo_defect =(EdgeVelo-MeanVelo_Profile)/v_star;
y_norm = height/BL_height;

%uncertainty in the velocity defect
norm_velo_defect_max =(EdgeVelo-MeanVelo_Profile-MeanVelo_Profile_uncertainty)/v_star;
norm_velo_defect_min =(EdgeVelo-MeanVelo_Profile+MeanVelo_Profile_uncertainty)/v_star;

end