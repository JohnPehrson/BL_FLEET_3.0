function [g1_bounds,g2_bounds] = widthbound_fy(run,ypix,gate1_location,gate2_location,...
    cutoff_height_switchover,g1_half_width,g2_half_width,lam_run_binary)
%This function will provide bounds for a double-gaussian fitting algorithm
%as a function of height in pixels. This will assume a constant set of
%bounds that then drop down to zero to work with boundary layers.

%initializing locational data
g1_bounds = zeros(ypix,3); %LB,x0,UB
g2_bounds = zeros(ypix,3);
g2_bound_extend = 8;

%do gate 1 as a line
for i = 1:ypix
    g1_bounds(i,:) = [gate1_location-g1_half_width,gate1_location,gate1_location+g1_half_width];
end

%do gate 2 as a line in the 'freestream'
g2frs = [1:cutoff_height_switchover(1),cutoff_height_switchover(3):ypix];
for i = g2frs
    g2_bounds(i,:) = [gate2_location-g2_half_width,gate2_location,gate2_location+g2_half_width+g2_bound_extend];
end

% figure;
% plot(g2_bounds,1:ypix);
% set(gca, 'YDir','reverse');

%get a function for decreasing velocity in g2bl
g2bl = [cutoff_height_switchover(1):ypix];
rows = g2bl-cutoff_height_switchover(2);

for temp = 1:1
above_wall = rows<0;
below_wall = rows>0;
rows = rows+above_wall.*1-below_wall.*1;
end

if lam_run_binary %true if more laminar
    colpix = gate1_location+(gate2_location-gate1_location).*(abs(rows)./max(abs(rows)));
else  %more turbulent runs
    colpix = (gate2_location-gate1_location).*log(abs(rows))./log(abs(max(abs(rows)))) +gate1_location;  %logorithmic, better for turbulent data
end
colpix(colpix<gate1_location) = gate1_location;

for i = 1:length(g2bl)
    g2_bounds(g2bl(i),2:3) = [colpix(i),gate2_location/2+colpix(i)/2+g2_half_width+g2_bound_extend];
end

for temp = 1:6
above_wall = rows<0;
below_wall = rows>0;
rows = rows+above_wall.*1-below_wall.*1;
end

if  lam_run_binary %true if more laminar
    colpix = gate1_location+(gate2_location-gate1_location).*(abs(rows)./max(abs(rows)));
else  %more turbulent runs
    colpix = (gate2_location-gate1_location).*log(abs(rows))./log(abs(max(abs(rows)))) +gate1_location;  %logorithmic, better for turbulent data
end
colpix(colpix<gate1_location) = gate1_location;

for i = 1:length(g2bl)
    g2_bounds(g2bl(i),1) = [colpix(i)-g2_half_width];
end

if lam_run_binary 
    for i = 1:length(g2bl)
        g2_bounds(g2bl(i),3) = [colpix(i)+g2_half_width+g2_bound_extend];
    end
end


figure;
plot(g2_bounds,1:ypix);
set(gca, 'YDir','reverse');

end