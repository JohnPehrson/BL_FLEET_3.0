function VelocityBoxPlots(red_velocity,emissionlocatingdata)
%This function plots a histogram of the instantaneous velocity at a variety
%of heights above the surface

rows = linspace(50,emissionlocatingdata(2),10);
    %make a matrix for the box plots
    box_mat = fliplr(red_velocity(rows,:)');
    boxplot_names = cell(size(box_mat,2),1);

    for i = 1:length(boxplot_names)
        fliprows = fliplr(rows);
        boxplot_names{i} = ['Velocity on row ',num2str(fliprows(i))];
    end
    
figure;
boxplot(box_mat(:,1:length(rows)),'Notch','on','Labels',boxplot_names(1:length(rows)))
ylabel('Streamwise Pixels');
xlabel('Velocity Distribution')

view([90 -90])

end