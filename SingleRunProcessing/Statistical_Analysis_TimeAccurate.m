function [P_vec,H_vec,mean_rms_centroids_current] = Statistical_Analysis_TimeAccurate(red_centroids,SNR,y_mm,gate1_location_bounds,gate2_location_bounds,...
    imageData_mean,mean_rms_centroids_prev)
%This function will provide statistical and visual analytics to show how
%well the data is fit, and the impact of filtering 'poor data'.

c1 = red_centroids(:,1:2:end);
c2 = red_centroids(:,2:2:end);

total_rows = size(c1,1);
nbins = 25;
all_rows = 1:length(y_mm);
wall_row = all_rows(y_mm==0);
use_rows = all_rows(~isnan(max(c1')));

mean_c = zeros(total_rows,2);
rms_c = zeros(total_rows,2);
for i = use_rows
    non_nan_c1 = c1(i,~isnan(c1(i,:)));
    non_nan_c2 = c2(i,~isnan(c1(i,:)));
    mean_c(i,:) = [mean(non_nan_c1),mean(non_nan_c2)];
    rms_c(i,:) =  [std(non_nan_c1),std(non_nan_c2)];
end

%% SNR statistics
SNR_mean = zeros(total_rows,1);
SNR_median = zeros(total_rows,1);
SNR_std = zeros(total_rows,1);

for i = use_rows
    row_nonNanSNR = SNR(i,~isnan(SNR(i,:)));
SNR_mean(i) = mean(row_nonNanSNR);
SNR_median(i) = median(row_nonNanSNR);
SNR_std(i) = std(row_nonNanSNR);
end

%histograms at 3 heights, showing fitting bounds
rows_hist = [50,150,250];
if use_rows(1)>1
    rows_hist(1) = use_rows(1);
end

    figure;
    subplot(3,2,1);
    histogram(c1(rows_hist(1),:),nbins);
        hold on;
    xline(gate1_location_bounds(rows_hist(1),1),'r','Linewidth',2);
    xline(gate1_location_bounds(rows_hist(1),3),'r','Linewidth',2);
        xlabel('X Pix Loc');
        ylabel('Probability');
        set(gca,'FontSize', 20);
        set(gca,'fontname','times')  % Set it to times
        title(['C1 row number ',num2str(rows_hist(1))]);
    subplot(3,2,3);
    histogram(c1(rows_hist(2),:),nbins);
        hold on;
    xline(gate1_location_bounds(rows_hist(2),1),'r','Linewidth',2);
    xline(gate1_location_bounds(rows_hist(2),3),'r','Linewidth',2);
        xlabel('X Pix Loc');
        ylabel('Probability');
        set(gca,'FontSize', 20);
        set(gca,'fontname','times')  % Set it to times
        title(['C1 row number ',num2str(rows_hist(2))]);
    subplot(3,2,5);
    histogram(c1(rows_hist(3),:),nbins);
        hold on;
    xline(gate1_location_bounds(rows_hist(3),1),'r','Linewidth',2);
    xline(gate1_location_bounds(rows_hist(3),3),'r','Linewidth',2);
        xlabel('X Pix Loc');
        ylabel('Probability');
        set(gca,'FontSize', 20);
        set(gca,'fontname','times')  % Set it to times
        title(['C1 row number ',num2str(rows_hist(3))]);
    subplot(3,2,2);
    histogram(c2(rows_hist(1),:),nbins);
        hold on;
    xline(gate2_location_bounds(rows_hist(1),1),'r','Linewidth',2);
    xline(gate2_location_bounds(rows_hist(1),3),'r','Linewidth',2);
        xlabel('X Pix Loc');
        ylabel('Probability');
        set(gca,'FontSize', 20);
        set(gca,'fontname','times')  % Set it to times
        title(['C2 row number ',num2str(rows_hist(1))]);
    subplot(3,2,4);
    histogram(c2(rows_hist(2),:),nbins);
        hold on;
    xline(gate2_location_bounds(rows_hist(2),1),'r','Linewidth',2);
    xline(gate2_location_bounds(rows_hist(2),3),'r','Linewidth',2);
        xlabel('X Pix Loc');
        ylabel('Probability');
        set(gca,'FontSize', 20);
        set(gca,'fontname','times')  % Set it to times
        title(['C2 row number ',num2str(rows_hist(2))]);
    subplot(3,2,6);
    histogram(c2(rows_hist(3),:),nbins);
        hold on;
    xline(gate2_location_bounds(rows_hist(3),1),'r','Linewidth',2);
    xline(gate2_location_bounds(rows_hist(3),3),'r','Linewidth',2);
        xlabel('X Pix Loc');
        ylabel('Probability');
        set(gca,'FontSize', 20);
        set(gca,'fontname','times')  % Set it to times
        title(['C2 row number ',num2str(rows_hist(3))]);

%box plot at 3 heights
rows = [50:25:325,347];
    %make a matrix for the box plots
    box_mat = [fliplr(c1(rows,:)'),fliplr(c2(rows,:)')];
    
    boxplot_names = cell(size(box_mat,2)/2,1);

    for i = 1:length(boxplot_names)
        fliprows = fliplr(rows);
        boxplot_names{i} = ['Row ',num2str(fliprows(i))];
    end
    ub_c1 = flipud(gate1_location_bounds(rows,3));
    lb_c1 = flipud(gate1_location_bounds(rows,1));
    ub_c2 = flipud(gate2_location_bounds(rows,3));
    lb_c2 = flipud(gate2_location_bounds(rows,1));

figure;
boxplot(box_mat(:,1:length(lb_c1)),'Notch','on','Labels',boxplot_names(1:length(lb_c1)))
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),'g','FaceAlpha',.5);
    end

hold on;
boxplot(box_mat(:,(length(lb_c1)+1):(length(lb_c2)+length(lb_c1))),'Notch','on','Labels',boxplot_names(1:length(lb_c1)))
    h = findobj(gca,'Tag','Box');
    for j=1:(length(h)/2)
        patch(get(h(j),'XData'),get(h(j),'YData'),'b','FaceAlpha',.5);
    end

line(1:length(ub_c1),ub_c1,'Color','red','LineStyle','--');
line(1:length(lb_c1),lb_c1,'Color','red','LineStyle','--');
line(1:length(ub_c2),ub_c2,'Color','red');
line(1:length(lb_c2),lb_c2,'Color','red');
ylabel('Streamwise Pixels')
ylim([min(lb_c1)-2,max(ub_c2)+2])
set(gca,'FontSize', 15);

leg_labels = ["\mu_1","","","","","","","","","","","","","\mu_2","","","","","","","","","","","","","Bounds for \mu_1","","Bounds for \mu_2",""];
legend(leg_labels);


view([90 -90])

%% QQ plot at a given height
row_qq = 200;
c1_dist = c1(row_qq,:);
c1_dist = c1_dist(~isnan(c1_dist));
c1_norm = normalize(c1_dist);
c2_dist = c2(row_qq,:);
c2_dist = c2_dist(~isnan(c2_dist));
c2_norm = normalize(c2_dist);
c1_bounds = gate1_location_bounds(row_qq,[1,3]);
c2_bounds = gate2_location_bounds(row_qq,[1,3]);
c1_bounds_norm = (c1_bounds-mean(c1_dist))./std(c1_dist);
c2_bounds_norm = (c2_bounds-mean(c2_dist))./std(c2_dist);

figure;
qqplot(c1_dist);
hold on;
qqplot(c2_dist);
legend('C1','C2');
title(['QQ Plot at row ',num2str(row_qq)])

figure;
qqplot(c1_norm);
hold on;
qqplot(c2_norm);
xline([c1_bounds_norm],'r','Linewidth',2)
xline([c2_bounds_norm],'b','Linewidth',2)
legend('','','C1 Distribution','','','C2 Distribution','C1 Bounds','','C2 Bounds','');
title(['Normalized QQ Plot at row ',num2str(row_qq)])

%normality tests using a One-sample Kolmogorov-Smirnov test
P_vec = zeros(total_rows,2);
H_vec = zeros(total_rows,2);
for i = use_rows
    c1_dist_row = c1(i,:);
    c1_dist_row = c1_dist_row(~isnan(c1_dist_row));
    c1_norm_row = normalize(c1_dist_row);
    c2_dist_row = c2(i,:);
    c2_dist_row = c2_dist_row(~isnan(c2_dist_row));
    c2_norm_row = normalize(c2_dist_row);

    [h1,p1] = kstest(c1_norm_row);
    [h2,p2] = kstest(c2_norm_row);
    P_vec(i,:) = [p1,p2];
    H_vec(i,:) = [h1,h2];
end

figure;
subplot(1,2,1);
    semilogx(P_vec(:,1),1:total_rows,'r','Linewidth',2);
    hold on;
    semilogx(P_vec(:,2),1:total_rows,'b','Linewidth',2);
    xline(0.025,'k')
    yline(wall_row,'m');
    grid on;
    xlabel('P value showing how normally distributed the data is');
    ylabel('Image Pixel Row');
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times
    legend('C1','C2','95% Cutoff','Wall Location')
    title(['Normality of Data for all Rows']);
    set(gca, 'YDir','reverse')
    xlim([10^(-10),5])
    ylim([1,total_rows])

subplot(1,2,2);
    plot(SNR_median,1:total_rows,'r','Linewidth',2);
    hold on;
    plot(SNR_median-SNR_std,1:total_rows,':r','Linewidth',2);
    plot(SNR_median+SNR_std,1:total_rows,':r','Linewidth',2);
    yline(wall_row,'m');
    grid on;
    xlabel('SNR');
    ylabel('Image Pixel Row');
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times
    legend('SNR','+/- 1 standard deviation','','Wall Location')
    title(['SNR']);
    set(gca, 'YDir','reverse')
    ylim([1,total_rows])


%Plot residual to the mean through time at one height
time = 1:length(c1_dist);
c1_resid = c1_dist-mean(c1_dist);
c2_resid = c2_dist-mean(c2_dist);

figure;
plot(time,c2_resid,'r','Linewidth',2)
hold on;
plot(time,c1_resid,'b','Linewidth',2)
grid on;
xlabel('Sequential Image Number');
ylabel('Residual to the Mean at that Image');
set(gca,'FontSize', 20);
set(gca,'fontname','times')  % Set it to times
legend('C1','C2')
title(['Residual to the mean through time, row number ',num2str(rows(3))]);

%Plot mean and rms of centroid location
    figure;
subplot(1,2,1);
    plot(mean_c(:,1),1:total_rows,'r','Linewidth',2);
    hold on;
    plot(mean_c(:,2),1:total_rows,'b','Linewidth',2);
    yline(wall_row,'m');
    grid on;
    xlabel('Mean Location of Centroids');
    ylabel('Image Pixel Row');
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times
    legend('C1','C2')
    set(gca, 'YDir','reverse')
    ylim([1,total_rows])

subplot(1,2,2);
    plot(rms_c(:,1),1:total_rows,'r','Linewidth',2);
    hold on;
    plot(rms_c(:,2),1:total_rows,'b','Linewidth',2);
    yline(wall_row,'m');
    grid on;
    xlabel('RMS of Centroids');
    ylabel('Image Pixel Row');
    set(gca,'FontSize', 20);
    set(gca,'fontname','times')  % Set it to times
    legend('C1','C2')
    set(gca, 'YDir','reverse')
    ylim([1,total_rows])

    %c2, nearwall histogram
                figure;
        histogram(c2(wall_row,:),40);
        hold on;
        xline(gate2_location_bounds(wall_row,1),'r','Linewidth',2);
        xline(gate2_location_bounds(wall_row,3),'r','Linewidth',2);
        xlabel('X Pix Loc');
        ylabel('Probability');
        set(gca,'FontSize', 20);
        set(gca,'fontname','times')  % Set it to times
        title(['C2 histogram at the wall']);


mean_rms_centroids_current = [mean_c,rms_c];

    if size(mean_rms_centroids_prev,2) >1
        mean_c_prev = mean_rms_centroids_prev(:,1:2);
        rms_c_prev = mean_rms_centroids_prev(:,3:4);

        mean_diff = mean_c-mean_c_prev;
        rms_diff = rms_c-rms_c_prev;

        %Plot change in mean and rms of centroids
            figure;
        subplot(1,2,1);
            plot(mean_diff(:,1),1:total_rows,'r','Linewidth',2);
            hold on;
            plot(mean_diff(:,2),1:total_rows,'b','Linewidth',2);
            yline(wall_row,'m');
            grid on;
            xlabel('Change in Mean Location of Centroids');
            ylabel('Image Pixel Row');
            set(gca,'FontSize', 20);
            set(gca,'fontname','times')  % Set it to times
            legend('C1','C2')
            set(gca, 'YDir','reverse')
            ylim([1,total_rows])
            xlim([-2,2]);
        
        subplot(1,2,2);
            plot(rms_diff(:,1),1:total_rows,'r','Linewidth',2);
            hold on;
            plot(rms_diff(:,2),1:total_rows,'b','Linewidth',2);
            yline(wall_row,'m');
            grid on;
            xlabel('Change inRMS of Centroids');
            ylabel('Image Pixel Row');
            set(gca,'FontSize', 20);
            set(gca,'fontname','times')  % Set it to times
            legend('C1','C2')
            set(gca, 'YDir','reverse')
            ylim([1,total_rows])
            xlim([-5,5]);

        %c2, uppermost, histogram
        figure;
        histogram(c2(rows_hist(1),:),40);
        hold on;
    xline(gate2_location_bounds(rows_hist(1),1),'r','Linewidth',2);
    xline(gate2_location_bounds(rows_hist(1),3),'r','Linewidth',2);
        xlabel('X Pix Loc');
        ylabel('Probability');
        set(gca,'FontSize', 20);
        set(gca,'fontname','times')  % Set it to times
        title(['C2 row number ',num2str(rows_hist(1))]);

    end

    
end