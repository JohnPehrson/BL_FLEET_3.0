function [uncertainty_mean_velo_wall_loc,height_unc] = Wall_Uncertainty_Propogator(resolution,wall_loc_unc,mean_velocity,heights)
%This function will propogate the uncertainty from the wall location into
%the mean velocity

            %get uncertainty in height in millimeters
            height_unc = wall_loc_unc.*(resolution(1)/(10^3));

            %move height up and down
                heights_dn = heights-height_unc;
                heights_up = heights+height_unc;
            %interpolate results if moved up or down
                velo_dn = interp1(heights_dn,mean_velocity,heights,'spline');
                velo_up = interp1(heights_up,mean_velocity,heights,'spline');
            %provide result
                uncertainty_mean_velo_wall_loc = (abs(velo_dn-mean_velocity)+abs(velo_up-mean_velocity))/2;

%             %plot
%                 figure;
%                 subplot(1,2,1);
%                 plot(mean_velocity,heights,'b','Linewidth',2);
%                 hold on;
%                 plot(mean_velocity,heights_dn);
%                 plot(mean_velocity,heights_up);
%                 plot(velo_dn,heights,':r','Linewidth',1);
%                 plot(velo_up,heights,':r','Linewidth',1);
%                     plot(mean_velocity-uncertainty_mean_velo_wall_loc,heights,':b','Linewidth',2);
%                     plot(mean_velocity+uncertainty_mean_velo_wall_loc,heights,':b','Linewidth',2);
%                 grid on;
%                 xlabel('Mean Velocity [m/s]')
%                 ylabel('Height above the surface [mm]')
%                 title('Uncertainty due to wall location Finding')
% 
%                 subplot(1,2,2);
%                 plot(uncertainty_mean_velo_wall_loc,heights,'b','Linewidth',2);
%                 xlim([0,50]);
%                 xlabel('Velocity Uncertainty [m/s]')
%                 title('Uncertainty due to wall location Finding')


end