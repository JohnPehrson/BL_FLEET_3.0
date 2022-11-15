function [resolution_height_uncertainty,wall_location_height_uncertainty] ...
    = Wall_Resolution_height_Uncertainty_Propogator(resolution,wall_loc_unc,heights)
%This function will propogate the uncertainty from the wall location into
%the mean velocity


    %% Propogating uncertainty in resolution into the height
        Pixel_vector = round(heights./(resolution(1)/(10^3)));
        res_unc = (resolution(2)/(10^3));
        resolution_height_uncertainty = Pixel_vector.*res_unc;
        
    %% Propogating the uncertianty in the wall location into the height
        wall_location_height_uncertainty = wall_loc_unc.*(resolution(1)/(10^3));


end