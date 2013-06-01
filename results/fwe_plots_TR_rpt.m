clear all
clc
close all

% This corresponds to 8 mm of smoothing in both cases
% SPM 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 mm
% RPT 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 mm
smoothing_plot = 5;

spm_datasets_voxel = zeros(7,8,3);
spm_datasets_cluster = zeros(7,8,3);

sums_rpt_voxel = zeros(7,8,3);
sums_rpt_cluster = zeros(7,8,3);

load rest_TRs.mat

files = 1:1484;
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
            
    combination = 'rpt_8_mm_1000_perm_AR4_cluster';
    
    % Voxel
    filename = ['../' combination '/max_test_values_rpt_' experiment_name];
    load(filename)
    filename = ['../' combination '/significant_voxels_rpt_' experiment_name];
    load(filename)
    filename = ['../' combination '/voxel_thresholds_rpt_' experiment_name];
    load(filename)
    
    smoothings = [1 3 5 7 9 11 13]; 
    for f = files
        if (abs(round(rest_TRs(f)) - rest_TRs(f)) < 0.05)
            for smoothing_ = 1:7
                smoothing = smoothings(smoothing_);
                spm_datasets_voxel(smoothing_,experiment,round(rest_TRs(f))) = spm_datasets_voxel(smoothing_,experiment,round(rest_TRs(f))) + 1;
                % The rpt script classifies datasets with a null threshold
                % as significant
                if voxel_thresholds_rpt(f,smoothing) > 0
                    sums_rpt_voxel(smoothing_,experiment,round(rest_TRs(f))) = sums_rpt_voxel(smoothing_,experiment,round(rest_TRs(f))) + significant_voxels_rpt(f,smoothing);                    
                end
            end
        end
    end
    
    % Cluster
    filename = ['../' combination '/max_cluster_sizes_rpt_' experiment_name];
    load(filename)
    filename = ['../' combination '/significant_clusters_rpt_' experiment_name];
    load(filename)
    filename = ['../' combination '/cluster_thresholds_rpt_' experiment_name];
    load(filename)
    
    smoothings = [1 3 5 7 9 11 13];
    for f = files
        if (abs(round(rest_TRs(f)) - rest_TRs(f)) < 0.05)
            for smoothing_ = 1:7
                smoothing = smoothings(smoothing_);
                spm_datasets_cluster(smoothing_,experiment,round(rest_TRs(f))) = spm_datasets_cluster(smoothing_,experiment,round(rest_TRs(f))) + 1;
                % The rpt script classifies datasets with a null threshold
                % as significant
                if cluster_thresholds_rpt(f,smoothing) > 0
                    sums_rpt_cluster(smoothing_,experiment,round(rest_TRs(f))) = sums_rpt_cluster(smoothing_,experiment,round(rest_TRs(f))) + significant_clusters_rpt(f,smoothing);                    
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
            voxel_FWEs(smoothing,experiment,TR) = sums_rpt_voxel(smoothing,experiment,TR) / spm_datasets_voxel(smoothing,experiment,TR);
            cluster_FWEs(smoothing,experiment,TR) = sums_rpt_cluster(smoothing,experiment,TR) / spm_datasets_cluster(smoothing,experiment,TR);
        end
    end
end




figure(1)
plot(5*ones(24,1),':k')
hold on
plot([1*ones(8,1); 3.5*ones(8,1); 2*ones(8,1)],'-.k','Color',[0.5 0.5 0.5])
hold on
plot(voxel_FWEs(5,:)*100,'-k')
hold on
plot([9*ones(8,1); 6.5*ones(8,1); 8*ones(8,1)],'-.k','Color',[0.5 0.5 0.5])
hold off
title('Voxel level inference, random permutation test','FontSize',15)
xlabel('TR = 1 s                   TR = 2 s                   TR = 3 s','FontSize',15)
ylabel('Familywise error rate (%)','FontSize',15)
legend('Truth','95% CI')
ylim([0 80])
set(gca,'FontSize',13)
L = [1 24];
NumTicks = 24;
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
set(gca,'XTickLabel',{'B1', 'B2', 'B3', 'B4', 'E1', 'E2', 'E3', 'E4' ,'B1', 'B2', 'B3', 'B4' ,'E1', 'E2', 'E3' ,'E4', 'B1' ,'B2', 'B3', 'B4' ,'E1' ,'E2', 'E3' ,'E4'})



print -dpng 'fwecomparison_rpt_voxel.png'
print -deps 'fwecomparison_rpt_voxel.eps'


figure(2)
plot(5*ones(24,1),':k')
hold on
plot([1*ones(8,1); 3.5*ones(8,1); 2*ones(8,1)],'-.k','Color',[0.5 0.5 0.5])
hold on
plot(cluster_FWEs(5,:)*100,'-k')
hold on
plot([9*ones(8,1); 6.5*ones(8,1); 8*ones(8,1)],'-.k','Color',[0.5 0.5 0.5])
hold off
title('Cluster level inference, random permutation test','FontSize',15)
xlabel('TR = 1 s                   TR = 2 s                   TR = 3 s','FontSize',15)
ylabel('Familywise error rate (%)','FontSize',15)
legend('Truth','95% CI')
ylim([0 80])
set(gca,'FontSize',13)
L = [1 24];
NumTicks = 24;
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
set(gca,'XTickLabel',{'B1', 'B2', 'B3', 'B4', 'E1', 'E2', 'E3', 'E4' ,'B1', 'B2', 'B3', 'B4' ,'E1', 'E2', 'E3' ,'E4', 'B1' ,'B2', 'B3', 'B4' ,'E1' ,'E2', 'E3' ,'E4'})


print -dpng 'fwecomparison_rpt_cluster.png'
print -deps 'fwecomparison_rpt_cluster.eps'


