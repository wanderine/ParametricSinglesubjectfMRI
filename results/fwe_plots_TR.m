clear all
clc
close all

% This corresponds to 8 mm of smoothing in both cases
% SPM 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 mm
% RPT 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 mm
smoothing_plot = 5;

spm_datasets_voxel = zeros(7,8,3);
spm_datasets_cluster = zeros(7,8,3);

sums_rft_voxel = zeros(7,8,3);
sums_rft_cluster = zeros(7,8,3);

load rest_TRs.mat


for experiment = 1:8
    
    if experiment == 1
        experiment_name = 'boxcar10';
    elseif experiment == 2
        experiment_name = 'boxcar15';
    elseif experiment == 3
        experiment_name = 'boxcar20';
    elseif experiment == 4
        experiment_name = 'boxcar30';
    elseif experiment == 5
        experiment_name = 'event1';
    elseif experiment == 6
        experiment_name = 'event2';
    elseif experiment == 7
        experiment_name = 'event3';
    elseif experiment == 8
        experiment_name = 'event4';
    end
    
   
    for file_range = 1:4
        
        if file_range == 1
            files = 1:375;
            files_string = '1_375';
        elseif file_range == 2
            files = 376:750;
            files_string = '376_750';
        elseif file_range == 3
            files = 751:1125;
            files_string = '751_1125';
        elseif file_range == 4
            files = 1126:1484;
            files_string = '1126_1484';
        end
        
        %combination = 'SPM_ufp0.05_scaling_no_motion_regressors_voxel_threshold_0.001';
        %combination = 'SPM_ufp0.05_no_scaling_motion_regressors_voxel_threshold_0.001';
        %combination = 'SPM_ufp0.05_no_scaling_no_motion_regressors_voxel_threshold_0.001';
        combination = 'SPM_ufp0.05_scaling_motion_regressors_voxel_threshold_0.001';
        
        % Voxel
        filename = ['../' combination '/max_test_values_spm_' experiment_name '_files' files_string];
        load(filename)
        filename = ['../' combination '/significant_voxels_spm_' experiment_name '_files' files_string];
        load(filename)
        filename = ['../' combination '/voxel_thresholds_spm_' experiment_name '_files' files_string];
        load(filename)
        
        smoothings = [1 3 5 7 9 11 13]; % Corresponds to 4 6 8 10 12 14 16 mm
        for f = files
            if (abs(round(rest_TRs(f)) - rest_TRs(f)) < 0.05)
                for smoothing_ = 1:7
                    smoothing = smoothings(smoothing_);
                    sums_rft_voxel(smoothing_,experiment,round(rest_TRs(f))) = sums_rft_voxel(smoothing_,experiment,round(rest_TRs(f))) + significant_voxels_spm(f,smoothing);                        
                    spm_datasets_voxel(smoothing_,experiment,round(rest_TRs(f))) = spm_datasets_voxel(smoothing_,experiment,round(rest_TRs(f))) + 1;
                end
            end
        end
        
        % Cluster
        filename = ['../' combination '/max_cluster_sizes_spm_' experiment_name '_files' files_string];
        load(filename)
        filename = ['../' combination '/significant_clusters_spm_' experiment_name '_files' files_string];
        load(filename)
        filename = ['../' combination '/cluster_thresholds_spm_' experiment_name '_files' files_string];
        load(filename)
        
        smoothings = [1 3 5 7 9 11 13]; % Corresponds to 4 6 8 10 12 14 16 mm
        for f = files
            if (abs(round(rest_TRs(f)) - rest_TRs(f)) < 0.05)
                for smoothing_ = 1:7
                    smoothing = smoothings(smoothing_);
                    sums_rft_cluster(smoothing_,experiment,round(rest_TRs(f))) = sums_rft_cluster(smoothing_,experiment,round(rest_TRs(f))) + significant_clusters_spm(f,smoothing);                        
                    spm_datasets_cluster(smoothing_,experiment,round(rest_TRs(f))) = spm_datasets_cluster(smoothing_,experiment,round(rest_TRs(f))) + 1;
                end                
            end
        end
        
    end
end

voxel_FWEs = zeros(7,8,3);
cluster_FWEs = zeros(7,8,3);
for experiment = 1:8
    for TR = 1:3
        for smoothing = 1:7
            voxel_FWEs(smoothing,experiment,TR) = sums_rft_voxel(smoothing,experiment,TR) / spm_datasets_voxel(smoothing,experiment,TR);
            cluster_FWEs(smoothing,experiment,TR) = sums_rft_cluster(smoothing,experiment,TR) / spm_datasets_cluster(smoothing,experiment,TR);
        end
    end
end


voxel_FWEs = reshape(voxel_FWEs,7,24);
cluster_FWEs = reshape(cluster_FWEs,7,24);

black = 0;
colour = 1;

if black == 1

    figure(1)
    plot(5*ones(24,1),':k')
    hold on
    plot([1*ones(8,1); 3.5*ones(8,1); 2*ones(8,1)],'-.k','Color',[0.5 0.5 0.5])
    hold on
    plot([9*ones(8,1); 6.5*ones(8,1); 8*ones(8,1)],'-.k','Color',[0.5 0.5 0.5])
    hold on 
    plot(voxel_FWEs(1,:)*100,'k')
    for smoothing = 2:7
        hold on
        plot(voxel_FWEs(smoothing,:)*100,'k')
    end    
    hold off
    %title('Voxel level inference, SPM8, no global normalization, no motion regressors','FontSize',15)
    title('Voxel level inference, SPM8, global normalization, motion regressors','FontSize',15)
    xlabel('TR = 1 s                   TR = 2 s                   TR = 3 s','FontSize',15)
    ylabel('Familywise error rate (%)','FontSize',15)
    legend('Truth','95% confidence interval')
    ylim([0 80])
    set(gca,'FontSize',13)
    NumTicks = 24;
    L = [1 24];
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    set(gca,'XTickLabel',{'B1', 'B2', 'B3', 'B4', 'E1', 'E2', 'E3', 'E4' ,'B1', 'B2', 'B3', 'B4' ,'E1', 'E2', 'E3' ,'E4', 'B1' ,'B2', 'B3', 'B4' ,'E1' ,'E2', 'E3' ,'E4'})

elseif colour == 1
    
    figure(1)
    plot(voxel_FWEs(1,:)*100,'b')
    hold on
    plot(voxel_FWEs(2,:)*100,'g')
    hold on
    plot(voxel_FWEs(3,:)*100,'r')
    hold on
    plot(voxel_FWEs(4,:)*100,'c')
    hold on
    plot(voxel_FWEs(5,:)*100,'m')
    hold on
    plot(voxel_FWEs(6,:)*100,'y')
    hold on
    plot(voxel_FWEs(7,:)*100,':k')
    hold on
    plot(5*ones(24,1),'k')
    hold on
    plot([1*ones(8,1); 3.5*ones(8,1); 2*ones(8,1)],'-.k','Color',[0.5 0.5 0.5])
    hold on
    plot([9*ones(8,1); 6.5*ones(8,1); 8*ones(8,1)],'-.k','Color',[0.5 0.5 0.5])
    hold off
    %title('Voxel level inference, SPM8, no global normalization, no motion regressors','FontSize',15)
    title('Voxel level inference, SPM8, global normalization, motion regressors','FontSize',15)
    xlabel('TR = 1 s                   TR = 2 s                   TR = 3 s','FontSize',15)
    ylabel('Familywise error rate (%)','FontSize',15)
    %legend('Truth','95% confidence interval')
    legend('4 mm','6 mm','8 mm','10 mm','12 mm','14 mm','16 mm','Truth','95% CI')
    ylim([0 80])
    set(gca,'FontSize',13)
    NumTicks = 24;
    L = [1 24];
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    set(gca,'XTickLabel',{'B1', 'B2', 'B3', 'B4', 'E1', 'E2', 'E3', 'E4' ,'B1', 'B2', 'B3', 'B4' ,'E1', 'E2', 'E3' ,'E4', 'B1' ,'B2', 'B3', 'B4' ,'E1' ,'E2', 'E3' ,'E4'})
    
end

%print -dpng 'fwecomparison_no_scaling_no_motion_regressors_voxel_black.png'
%print -deps 'fwecomparison_no_scaling_no_motion_regressors_voxel_black.eps'

%print -dpng 'fwecomparison_no_scaling_no_motion_regressors_voxel_colour.png'
%print -depsc 'fwecomparison_no_scaling_no_motion_regressors_voxel_colour.eps'

%print -dpng 'fwecomparison_scaling_motion_regressors_voxel_black.png'
%print -deps 'fwecomparison_scaling_motion_regressors_voxel_black.eps'

%print -dpng 'fwecomparison_scaling_motion_regressors_voxel_colour.png'
%print -depsc 'fwecomparison_scaling_motion_regressors_voxel_colour.eps'

if black == 1
    
    figure(2)
    plot(5*ones(24,1),':k')
    hold on
    plot([1*ones(8,1); 3.5*ones(8,1); 2*ones(8,1)],'-.k','Color',[0.5 0.5 0.5])
    hold on
    plot([9*ones(8,1); 6.5*ones(8,1); 8*ones(8,1)],'-.k','Color',[0.5 0.5 0.5])
    hold on
    plot(cluster_FWEs(1,:)*100,'k')
    for smoothing = 2:7
        hold on
        plot(cluster_FWEs(smoothing,:)*100,'k')
    end    
    hold off
    %title('Cluster level inference, SPM8, no global normalization, no motion regressors','FontSize',15)
    title('Cluster level inference, SPM8, global normalization, motion regressors','FontSize',15)
    xlabel('TR = 1 s                   TR = 2 s                   TR = 3 s','FontSize',15)
    ylabel('Familywise error rate (%)','FontSize',15)
    legend('Truth','95% confidence interval')
    ylim([0 80])
    set(gca,'FontSize',13)
    NumTicks = 24;
    L = [1 24];
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    set(gca,'XTickLabel',{'B1', 'B2', 'B3', 'B4', 'E1', 'E2', 'E3', 'E4' ,'B1', 'B2', 'B3', 'B4' ,'E1', 'E2', 'E3' ,'E4', 'B1' ,'B2', 'B3', 'B4' ,'E1' ,'E2', 'E3' ,'E4'})
    
    
elseif colour == 1
    
    figure(2)
    plot(cluster_FWEs(1,:)*100,'b')
    hold on
    plot(cluster_FWEs(2,:)*100,'g')
    hold on
    plot(cluster_FWEs(3,:)*100,'r')
    hold on
    plot(cluster_FWEs(4,:)*100,'c')
    hold on
    plot(cluster_FWEs(5,:)*100,'m')
    hold on
    plot(cluster_FWEs(6,:)*100,'y')
    hold on
    plot(cluster_FWEs(7,:)*100,':k')
    hold on
    plot(5*ones(24,1),'k')
    hold on
    plot([1*ones(8,1); 3.5*ones(8,1); 2*ones(8,1)],'-.k','Color',[0.5 0.5 0.5])
    hold on
    plot([9*ones(8,1); 6.5*ones(8,1); 8*ones(8,1)],'-.k','Color',[0.5 0.5 0.5])
    hold off
    %title('Cluster level inference, SPM8, no global normalization, no motion regressors','FontSize',15)
    title('Cluster level inference, SPM8, global normalization, motion regressors','FontSize',15)
    xlabel('TR = 1 s                   TR = 2 s                   TR = 3 s','FontSize',15)
    ylabel('Familywise error rate (%)','FontSize',15)
    %legend('Truth','95% confidence interval')
    legend('4 mm','6 mm','8 mm','10 mm','12 mm','14 mm','16 mm','Truth','95% CI')
    ylim([0 80])
    set(gca,'FontSize',13)
    NumTicks = 24;
    L = [1 24];
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    set(gca,'XTickLabel',{'B1', 'B2', 'B3', 'B4', 'E1', 'E2', 'E3', 'E4' ,'B1', 'B2', 'B3', 'B4' ,'E1', 'E2', 'E3' ,'E4', 'B1' ,'B2', 'B3', 'B4' ,'E1' ,'E2', 'E3' ,'E4'})
    
end


%print -dpng 'fwecomparison_no_scaling_no_motion_regressors_cluster_black.png'
%print -deps 'fwecomparison_no_scaling_no_motion_regressors_cluster_black.eps'

%print -dpng 'fwecomparison_no_scaling_no_motion_regressors_cluster_colour.png'
%print -depsc 'fwecomparison_no_scaling_no_motion_regressors_cluster_colour.eps'

%print -dpng 'fwecomparison_scaling_motion_regressors_cluster_black.png'
%print -deps 'fwecomparison_scaling_motion_regressors_cluster_black.eps'

%print -dpng 'fwecomparison_scaling_motion_regressors_cluster_colour.png'
%print -depsc 'fwecomparison_scaling_motion_regressors_cluster_colour.eps'


