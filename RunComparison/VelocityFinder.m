function [velocity,random_velocity_error,systematic_velocity_error] = VelocityFinder(centroids,centroid_error,resolution,gates,delays)
%This function finds the velocity of the flow based on the centroids.

        resolution(1) = resolution(1)/(10^6);
        t = [ delays(1)+gates(1)/2 ,delays(1)+gates(1)+delays(2)+gates(2)/2];
        dt = (t(2)-t(1))*10^(-9); %in seconds
        dt_err = 0; %no jitter between gates in a burst
       
        dx_pix = centroids(:,2)-centroids(:,1);
        dx = dx_pix.*resolution(1); %displacement in meters
        velocity = dx./dt;
        dx_pix_err = sqrt(centroid_error(:,1).^2+centroid_error(:,2).^2);
        
        reserr = resolution(2)/(10^5);
        
        random_du2 = (((resolution(1)/(dt)).^2).*(dx_pix_err).^2)+ (((dx_pix.*resolution(1)./(dt.^2)).^2).*(dt_err).^2);
        random_velocity_error = sqrt(random_du2);
        
        systematic_du2 = (((dx_pix./(dt)).^2).*(reserr)^2);
        systematic_velocity_error = sqrt(systematic_du2);
        
end

