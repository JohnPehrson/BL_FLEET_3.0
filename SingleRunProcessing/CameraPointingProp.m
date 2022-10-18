function [row_mm_out,row_mm_uncertainty_out,mean_velocity_s_out,rms_velocity_s_out] = CameraPointingProp(mean_velocity,...
    rms_velocity,mean_velocity_s,rms_velocity_s,inclination_angle,inclination_angle_unc,...
    emissionlocatingdata,pixel_um_resolution,row_mm,zero_height_ref_unc)
%This function changes the inclination and span of the image angle. 

%% Wall location uncertainty
loc_unc = (zero_height_ref_unc).*pixel_um_resolution(1)./1000;

%% Inclination angle
row_mm_out = row_mm./cosd(inclination_angle);
row_mm_uncertainty_int = abs(row_mm_out.*(sind(inclination_angle)./cosd(inclination_angle))).*deg2rad(inclination_angle_unc);

%inclination uncertainty + general camera scale uncertainty + wall location
%uncertainty
row_mm_uncertainty_out = sqrt(row_mm_uncertainty_int.^2+(pixel_um_resolution(2).*row_mm_out./1000).^2+loc_unc^2);

%% Span angle
%span angle ideal is at 0 degrees, so mean and rms won't change, just
%change systematic errors
spanangle_uncertainty = 0.5; %degree
mean_velocity_s_out = sqrt((cosd(0).^2).*(mean_velocity_s.^2) +(mean_velocity.^2).*(sind(spanangle_uncertainty).^2));
rms_velocity_s_out = sqrt((cosd(0).^2).*(rms_velocity_s.^2) +(rms_velocity.^2).*(sind(spanangle_uncertainty).^2));
end