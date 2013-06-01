%  
% max_test_values_rpt = zeros(1484,16);
% save max_test_values_rpt_boxcar10.mat max_test_values_rpt
% save max_test_values_rpt_boxcar15.mat max_test_values_rpt
% save max_test_values_rpt_boxcar20.mat max_test_values_rpt
% save max_test_values_rpt_boxcar30.mat max_test_values_rpt
% save max_test_values_rpt_event1.mat max_test_values_rpt
% save max_test_values_rpt_event2.mat max_test_values_rpt
% save max_test_values_rpt_event3.mat max_test_values_rpt
% save max_test_values_rpt_event4.mat max_test_values_rpt
% 
% max_cluster_sizes_rpt = zeros(1484,16);
% save max_cluster_sizes_rpt_boxcar10.mat max_cluster_sizes_rpt
% save max_cluster_sizes_rpt_boxcar15.mat max_cluster_sizes_rpt
% save max_cluster_sizes_rpt_boxcar20.mat max_cluster_sizes_rpt
% save max_cluster_sizes_rpt_boxcar30.mat max_cluster_sizes_rpt
% save max_cluster_sizes_rpt_event1.mat max_cluster_sizes_rpt
% save max_cluster_sizes_rpt_event2.mat max_cluster_sizes_rpt
% save max_cluster_sizes_rpt_event3.mat max_cluster_sizes_rpt
% save max_cluster_sizes_rpt_event4.mat max_cluster_sizes_rpt
% 
% 
% significant_voxels_rpt = zeros(1484,16);
% save significant_voxels_rpt_boxcar10.mat significant_voxels_rpt
% save significant_voxels_rpt_boxcar15.mat significant_voxels_rpt
% save significant_voxels_rpt_boxcar20.mat significant_voxels_rpt
% save significant_voxels_rpt_boxcar30.mat significant_voxels_rpt
% save significant_voxels_rpt_event1.mat significant_voxels_rpt
% save significant_voxels_rpt_event2.mat significant_voxels_rpt
% save significant_voxels_rpt_event3.mat significant_voxels_rpt
% save significant_voxels_rpt_event4.mat significant_voxels_rpt
% 
% significant_clusters_rpt = zeros(1484,16);
% save significant_clusters_rpt_boxcar10.mat significant_clusters_rpt
% save significant_clusters_rpt_boxcar15.mat significant_clusters_rpt
% save significant_clusters_rpt_boxcar20.mat significant_clusters_rpt
% save significant_clusters_rpt_boxcar30.mat significant_clusters_rpt
% save significant_clusters_rpt_event1.mat significant_clusters_rpt
% save significant_clusters_rpt_event2.mat significant_clusters_rpt
% save significant_clusters_rpt_event3.mat significant_clusters_rpt
% save significant_clusters_rpt_event4.mat significant_clusters_rpt
% 
% 
% 
% voxel_thresholds_rpt = zeros(1484,16);
% save voxel_thresholds_rpt_boxcar10.mat voxel_thresholds_rpt
% save voxel_thresholds_rpt_boxcar15.mat voxel_thresholds_rpt
% save voxel_thresholds_rpt_boxcar20.mat voxel_thresholds_rpt
% save voxel_thresholds_rpt_boxcar30.mat voxel_thresholds_rpt
% save voxel_thresholds_rpt_event1.mat voxel_thresholds_rpt
% save voxel_thresholds_rpt_event2.mat voxel_thresholds_rpt
% save voxel_thresholds_rpt_event3.mat voxel_thresholds_rpt
% save voxel_thresholds_rpt_event4.mat voxel_thresholds_rpt
% 
% cluster_thresholds_rpt = zeros(1484,16);
% save cluster_thresholds_rpt_boxcar10.mat cluster_thresholds_rpt
% save cluster_thresholds_rpt_boxcar15.mat cluster_thresholds_rpt
% save cluster_thresholds_rpt_boxcar20.mat cluster_thresholds_rpt
% save cluster_thresholds_rpt_boxcar30.mat cluster_thresholds_rpt
% save cluster_thresholds_rpt_event1.mat cluster_thresholds_rpt
% save cluster_thresholds_rpt_event2.mat cluster_thresholds_rpt
% save cluster_thresholds_rpt_event3.mat cluster_thresholds_rpt
% save cluster_thresholds_rpt_event4.mat cluster_thresholds_rpt


clear all
close all
clc

load rand_events1.mat
onsets1 = onsets;
durations1 = durations;

load rand_events2.mat
onsets2 = onsets;
durations2 = durations;

addpath('/home/andek/Testsaker/spm8');
addpath('/home/andek/Testsaker/nifti_matlab');
number_of_permutations = 1000;
alpha_FWHM = 8.01;
T_FWHM = 8.01;
smoothing_filter_size = 9;
plott = 1;

smoothing_filters_x = zeros(smoothing_filter_size,16);
smoothing_filters_y = zeros(smoothing_filter_size,16);
smoothing_filters_z = zeros(smoothing_filter_size,16);

for file_number = 1:1484

    file_number
    
    % Skip datasets with empty brain masks
    if (file_number == 905)
        file_number = 906;
    elseif (file_number == 1310)
        file_number = 1311;
    end
    
    % Load data
    filename = ['/flush/andek/rest_fMRI/func/rest_' num2str(file_number) '.nii'];
    V = spm_vol(filename);
    
    % Get number of time points
    st = size(V,1);
    st
    
    % Get repetition time
    a = V.private;
    TR = a.timing.tspace;
    
    % Get voxel size
    voxel_size_x = abs(a.mat(1,1));
    voxel_size_y = abs(a.mat(2,2));
    voxel_size_z = abs(a.mat(3,3));
    
    % Smoothing for thresholding for segmentation
    sigma_x = T_FWHM / 2.354 / voxel_size_x;
    sigma_y = T_FWHM / 2.354 / voxel_size_y;
    sigma_z = T_FWHM / 2.354 / voxel_size_z;
    
    smoothing_filter_x = fspecial('gaussian',[9 9],sigma_x);
    smoothing_filter_x = smoothing_filter_x(:,5);
    smoothing_filter_x = smoothing_filter_x / sum(abs(smoothing_filter_x(:)));
    
    smoothing_filter_y = fspecial('gaussian',[9 9],sigma_y);
    smoothing_filter_y = smoothing_filter_y(:,5);
    smoothing_filter_y = smoothing_filter_y / sum(abs(smoothing_filter_y(:)));
    
    smoothing_filter_z = fspecial('gaussian',[9 9],sigma_z);
    smoothing_filter_z = smoothing_filter_z(:,5);
    smoothing_filter_z = smoothing_filter_z / sum(abs(smoothing_filter_z(:)));
    
    x_smoothing_filter = zeros(1,smoothing_filter_size,1);
    x_smoothing_filter(1,:,1) = smoothing_filter_x;
    
    y_smoothing_filter = zeros(smoothing_filter_size,1,1);
    y_smoothing_filter(:,1,1) = smoothing_filter_y;
    
    z_smoothing_filter = zeros(1,1,smoothing_filter_size);
    z_smoothing_filter(1,1,:) = smoothing_filter_z;
    
    % Get data
    fMRI_volumes = load_untouch_nii(filename);
    fMRI_volumes = fMRI_volumes.img;
    fMRI_volumes = double(fMRI_volumes);    
    
    %size(fMRI_volumes)
    
    if ( (file_number >= 924) && (file_number <= 969) )
        fMRI_volumes = permute(fMRI_volumes,[2 3 1 4]);
    end
        
    %size(fMRI_volumes)
    
    volume = fMRI_volumes(:,:,:,1);
    smoothed_volume = conv3(volume + eps,x_smoothing_filter);
    smoothed_volume = conv3(smoothed_volume,y_smoothing_filter);
    smoothed_volume = conv3(smoothed_volume,z_smoothing_filter);
    brain_voxels = double(smoothed_volume > 0.8*mean(smoothed_volume(:)));
    certainty = brain_voxels;
    
    % Alpha smoothing (smoothing of AR estimates)
    
    sigma_x = alpha_FWHM / 2.354 / voxel_size_x;
    sigma_y = alpha_FWHM / 2.354 / voxel_size_y;
    sigma_z = alpha_FWHM / 2.354 / voxel_size_z;
    
    alpha_smoothing_filter_x = fspecial('gaussian',[9 9],sigma_x);
    alpha_smoothing_filter_x = alpha_smoothing_filter_x(:,5);
    alpha_smoothing_filter_x = alpha_smoothing_filter_x / sum(abs(alpha_smoothing_filter_x(:)));
    
    alpha_smoothing_filter_y = fspecial('gaussian',[9 9],sigma_y);
    alpha_smoothing_filter_y = alpha_smoothing_filter_y(:,5);
    alpha_smoothing_filter_y = alpha_smoothing_filter_y / sum(abs(alpha_smoothing_filter_y(:)));
    
    alpha_smoothing_filter_z = fspecial('gaussian',[9 9],sigma_z);
    alpha_smoothing_filter_z = alpha_smoothing_filter_z(:,5);
    alpha_smoothing_filter_z = alpha_smoothing_filter_z / sum(abs(alpha_smoothing_filter_z(:)));                                
        
    x_alpha_smoothing_filter = zeros(1,smoothing_filter_size,1);
    x_alpha_smoothing_filter(1,:,1) = alpha_smoothing_filter_x;
    
    y_alpha_smoothing_filter = zeros(smoothing_filter_size,1,1);    
    y_alpha_smoothing_filter(:,1,1) = alpha_smoothing_filter_y;
    
    z_alpha_smoothing_filter = zeros(1,1,smoothing_filter_size);
    z_alpha_smoothing_filter(1,1,:) = alpha_smoothing_filter_z;
    
    smoothed_alpha_certainty = conv3(brain_voxels + eps,x_alpha_smoothing_filter);
    smoothed_alpha_certainty = conv3(smoothed_alpha_certainty,y_alpha_smoothing_filter);
    smoothed_alpha_certainty = conv3(smoothed_alpha_certainty,z_alpha_smoothing_filter);
        
    [sy sx sz st] = size(fMRI_volumes);       
    
    % Calculate file size for 32 bit floats, in MB
    file_size = sx * sy * sz * st * 4 / 1024 / 1024;
    
    % Only use the GPU with 3 GB of memory
    if (file_size > 200)
        number_of_GPUs = 1;
    % Use all the GPUs
    else
        number_of_GPUs = 2;
    end
    
    % Detrending, cubic
    X_Detrend = zeros(st,4);
    X_Detrend(:,1) = ones(st,1);
    X_Detrend(:,2) = -(st-1)/2:(st-1)/2;
    X_Detrend(:,3) = X_Detrend(:,2).^2;
    X_Detrend(:,4) = X_Detrend(:,2).^3;
    
    X_Detrend(:,1) = X_Detrend(:,1) / norm(X_Detrend(:,1));
    X_Detrend(:,2) = X_Detrend(:,2) / norm(X_Detrend(:,2));
    X_Detrend(:,3) = X_Detrend(:,3) / norm(X_Detrend(:,3));
    X_Detrend(:,4) = X_Detrend(:,4) / norm(X_Detrend(:,4));
    
    xtxxt_Detrend = inv(X_Detrend'*X_Detrend)*X_Detrend';
    
    % Motion compensation
    load filters.mat
    motion_compensation_filter_size = 7;
    number_of_iterations = 5;
    
    % Smoothing
    smoothed_certainties = zeros(sy,sx,sz,16);
    for smoothing = 0:15;
        
        FWHM = smoothing + 0.01;
        sigma_x = FWHM / 2.354 / voxel_size_x;
        sigma_y = FWHM / 2.354 / voxel_size_y;
        sigma_z = FWHM / 2.354 / voxel_size_z;
        
        smoothing_filter_x = fspecial('gaussian',[9 9],sigma_x);
        smoothing_filter_x = smoothing_filter_x(:,5);
        smoothing_filter_x = smoothing_filter_x / sum(abs(smoothing_filter_x(:)));
        smoothing_filters_x(:,smoothing+1) = smoothing_filter_x;
        
        smoothing_filter_y = fspecial('gaussian',[9 9],sigma_y);
        smoothing_filter_y = smoothing_filter_y(:,5);
        smoothing_filter_y = smoothing_filter_y / sum(abs(smoothing_filter_y(:)));
        smoothing_filters_y(:,smoothing+1) = smoothing_filter_y;
        
        smoothing_filter_z = fspecial('gaussian',[9 9],sigma_z);
        smoothing_filter_z = smoothing_filter_z(:,5);
        smoothing_filter_z = smoothing_filter_z / sum(abs(smoothing_filter_z(:)));
        smoothing_filters_z(:,smoothing+1) = smoothing_filter_z;
        
        x_smoothing_filter = zeros(1,smoothing_filter_size,1);
        x_smoothing_filter(1,:,1) = smoothing_filter_x;
    
        y_smoothing_filter = zeros(smoothing_filter_size,1,1);
        y_smoothing_filter(:,1,1) = smoothing_filter_y;
    
        z_smoothing_filter = zeros(1,1,smoothing_filter_size);
        z_smoothing_filter(1,1,:) = smoothing_filter_z;
        
        smoothed_certainty = conv3(brain_voxels + eps,x_smoothing_filter);
        smoothed_certainty = conv3(smoothed_certainty,y_smoothing_filter);
        smoothed_certainty = conv3(smoothed_certainty,z_smoothing_filter);
        smoothed_certainties(:,:,:,smoothing+1) = smoothed_certainty;
        
    end
    
    
    % GLM        
    
    % 10 s boxcar
    X_GLM = zeros(st*10000,2);
    X_GLM_1 = zeros(st,2);    
    activity_length = 10;
    rest_length = 10;
    t = round((rest_length+1)/TR*10000); % start with rest
    activity = round((activity_length+1)/TR*10000);
    rest = round((rest_length+1)/TR*10000);
    while t < st*10000
        X_GLM(t:t+activity,1) = 1;        
        t = t + activity + rest;
    end
    X_GLM = X_GLM(1:st*10000,:);
    X_GLM_1(:,1) = decimate(X_GLM(:,1),10000);
    
    
    % 15 s boxcar
    X_GLM = zeros(st*10000,2);
    X_GLM_2 = zeros(st,2);    
    activity_length = 15;
    rest_length = 15;
    t = round((rest_length+1)/TR*10000); % start with rest
    activity = round((activity_length+1)/TR*10000);
    rest = round((rest_length+1)/TR*10000);
    while t < st*10000
        X_GLM(t:t+activity,1) = 1;        
        t = t + activity + rest;
    end
    X_GLM = X_GLM(1:st*10000,:);
    X_GLM_2(:,1) = decimate(X_GLM(:,1),10000);
    
    
    
    
    % 20 s boxcar
    X_GLM = zeros(st*10000,2);
    X_GLM_3 = zeros(st,2);    
    activity_length = 20;
    rest_length = 20;
    t = round((rest_length+1)/TR*10000); % start with rest
    activity = round((activity_length+1)/TR*10000);
    rest = round((rest_length+1)/TR*10000);
    while t < st*10000
        X_GLM(t:t+activity,1) = 1;        
        t = t + activity + rest;
    end
    X_GLM = X_GLM(1:st*10000,:);
    X_GLM_3(:,1) = decimate(X_GLM(:,1),10000);
    
    
    
    
    % 30 s boxcar
    X_GLM = zeros(st*10000,2);
    X_GLM_4 = zeros(st,2);    
    activity_length = 30;
    rest_length = 30;
    t = round((rest_length+1)/TR*10000); % start with rest
    activity = round((activity_length+1)/TR*10000);
    rest = round((rest_length+1)/TR*10000);
    while t < st*10000
        X_GLM(t:t+activity,1) = 1;        
        t = t + activity + rest;
    end
    X_GLM = X_GLM(1:st*10000,:);
    X_GLM_4(:,1) = decimate(X_GLM(:,1),10000);
    
    
    
    
    % Event 1, 2 s activity, 6 s rest    
    X_GLM = zeros(st*10000,2);
    X_GLM_5 = zeros(st,2);    
    activity_length = 2;
    rest_length = 6;
    t = round((rest_length+1)/TR*10000);
    activity = round((activity_length+1)/TR*10000);
    rest = round((rest_length+1)/TR*10000);
    while t < st*10000
        X_GLM(t:t+activity,1) = 1;        
        t = t + activity + rest;
    end
    X_GLM = X_GLM(1:st*10000,:);
    X_GLM_5(:,1) = decimate(X_GLM(:,1),10000);
    
    
    
    % Event 2, 4 s activity, 8 s rest    
    X_GLM = zeros(st*10000,2);
    X_GLM_6 = zeros(st,2);    
    activity_length = 4;
    rest_length = 8;
    t = round((rest_length+1)/TR*10000);    
    activity = round((activity_length+1)/TR*10000);
    rest = round((rest_length+1)/TR*10000);
    while t < st*10000
        X_GLM(t:t+activity,1) = 1;        
        t = t + activity + rest;
    end
    X_GLM = X_GLM(1:st*10000,:);
    X_GLM_6(:,1) = decimate(X_GLM(:,1),10000);
    
    
    
    
    
    % Calculate how many onsets and durations to use, 
    % from the TR and the number of scans
    stt = st * TR;
    last = 1;
    tot = 0;
    while 1
        tot = onsets1(last) + durations1(last);
        
        if tot > stt
            last = last - 1;
            tot = onsets1(last) + durations1(last);
            break
        end
        last = last + 1;        
    end
    
    % Event 3, rand 1
    X_GLM = zeros(st*10000,2);
    X_GLM_7 = zeros(st,2);
    activity = round((durations1(1)+1)/TR*10000);    
    t = round((onsets1(1)+1)/TR*10000);    
    for dur = 2:last
        X_GLM(t:t+activity,1) = 1;                
        t = round((onsets1(dur)+1)/TR*10000);    
        activity = round((durations1(dur)+1)/TR*10000);       
    end
    X_GLM = X_GLM(1:st*10000,:);
    X_GLM_7(:,1) = decimate(X_GLM(:,1),10000);
    
    
    
    
    % Calculate how many onsets and durations to use, 
    % from the TR and the number of scans
    stt = st * TR;
    last = 1;
    tot = 0;
    while 1
        tot = onsets2(last) + durations2(last);
        
        if tot > stt
            last = last - 1;
            tot = onsets2(last) + durations2(last);
            break
        end
        last = last + 1;        
    end
        
    % Event 4, rand 2
    X_GLM = zeros(st*10000,2);
    X_GLM_8 = zeros(st,2);
    activity = round((durations2(1)+1)/TR*10000);    
    t = round((onsets2(1)+1)/TR*10000);    
    for dur = 2:last
        X_GLM(t:t+activity,1) = 1;                
        t = round((onsets2(dur)+1)/TR*10000);    
        activity = round((durations2(dur)+1)/TR*10000);        
    end
    X_GLM = X_GLM(1:st*10000,:);
    X_GLM_8(:,1) = decimate(X_GLM(:,1),10000);
    
    
    

    
    
    % Convolve with HRF and temporal derivative
    
    hrf = spm_hrf(2);
    dhrf = diff(hrf);
    
    % 1
    GLM1 = conv(X_GLM_1(:,1),hrf);
    GLM1 = GLM1(1:st);
    
    GLM2 = conv(X_GLM_1(:,1),dhrf);
    GLM2 = GLM2(1:st);
    
    X_GLM_1(:,1) = GLM1 - mean(GLM1);
    X_GLM_1(:,2) = GLM2 - mean(GLM2);
    
    % 2
    GLM1 = conv(X_GLM_2(:,1),hrf);
    GLM1 = GLM1(1:st);
    
    GLM2 = conv(X_GLM_2(:,1),dhrf);
    GLM2 = GLM2(1:st);
    
    X_GLM_2(:,1) = GLM1 - mean(GLM1);
    X_GLM_2(:,2) = GLM2 - mean(GLM2);
    
    % 3
    GLM1 = conv(X_GLM_3(:,1),hrf);
    GLM1 = GLM1(1:st);
    
    GLM2 = conv(X_GLM_3(:,1),dhrf);
    GLM2 = GLM2(1:st);
    
    X_GLM_3(:,1) = GLM1 - mean(GLM1);
    X_GLM_3(:,2) = GLM2 - mean(GLM2);
    
    % 4
    GLM1 = conv(X_GLM_4(:,1),hrf);
    GLM1 = GLM1(1:st);
    
    GLM2 = conv(X_GLM_4(:,1),dhrf);
    GLM2 = GLM2(1:st);
    
    X_GLM_4(:,1) = GLM1 - mean(GLM1);
    X_GLM_4(:,2) = GLM2 - mean(GLM2);
    
    % 5
    GLM1 = conv(X_GLM_5(:,1),hrf);
    GLM1 = GLM1(1:st);
    
    GLM2 = conv(X_GLM_5(:,1),dhrf);
    GLM2 = GLM2(1:st);
    
    X_GLM_5(:,1) = GLM1 - mean(GLM1);
    X_GLM_5(:,2) = GLM2 - mean(GLM2);
    
    % 6
    GLM1 = conv(X_GLM_6(:,1),hrf);
    GLM1 = GLM1(1:st);
    
    GLM2 = conv(X_GLM_6(:,1),dhrf);
    GLM2 = GLM2(1:st);
    
    X_GLM_6(:,1) = GLM1 - mean(GLM1);
    X_GLM_6(:,2) = GLM2 - mean(GLM2);
    
    % 7
    GLM1 = conv(X_GLM_7(:,1),hrf);
    GLM1 = GLM1(1:st);
    
    GLM2 = conv(X_GLM_7(:,1),dhrf);
    GLM2 = GLM2(1:st);
    
    X_GLM_7(:,1) = GLM1 - mean(GLM1);
    X_GLM_7(:,2) = GLM2 - mean(GLM2);
    
    % 8
    GLM1 = conv(X_GLM_8(:,1),hrf);
    GLM1 = GLM1(1:st);
    
    GLM2 = conv(X_GLM_8(:,1),dhrf);
    GLM2 = GLM2(1:st);
    
    X_GLM_8(:,1) = GLM1 - mean(GLM1);
    X_GLM_8(:,2) = GLM2 - mean(GLM2);
    
    
    % Normalize
    X_GLM_1(:,1) = X_GLM_1(:,1)/norm(X_GLM_1(:,1));
    X_GLM_1(:,2) = X_GLM_1(:,2)/norm(X_GLM_1(:,2));
    
    X_GLM_2(:,1) = X_GLM_2(:,1)/norm(X_GLM_2(:,1));
    X_GLM_2(:,2) = X_GLM_2(:,2)/norm(X_GLM_2(:,2));
    
    X_GLM_3(:,1) = X_GLM_3(:,1)/norm(X_GLM_3(:,1));
    X_GLM_3(:,2) = X_GLM_3(:,2)/norm(X_GLM_3(:,2));
    
    X_GLM_4(:,1) = X_GLM_4(:,1)/norm(X_GLM_4(:,1));
    X_GLM_4(:,2) = X_GLM_4(:,2)/norm(X_GLM_4(:,2));
    
    X_GLM_5(:,1) = X_GLM_5(:,1)/norm(X_GLM_5(:,1));
    X_GLM_5(:,2) = X_GLM_5(:,2)/norm(X_GLM_5(:,2));
    
    X_GLM_6(:,1) = X_GLM_6(:,1)/norm(X_GLM_6(:,1));
    X_GLM_6(:,2) = X_GLM_6(:,2)/norm(X_GLM_6(:,2));
    
    X_GLM_7(:,1) = X_GLM_7(:,1)/norm(X_GLM_7(:,1));
    X_GLM_7(:,2) = X_GLM_7(:,2)/norm(X_GLM_7(:,2));
    
    X_GLM_8(:,1) = X_GLM_8(:,1)/norm(X_GLM_8(:,1));
    X_GLM_8(:,2) = X_GLM_8(:,2)/norm(X_GLM_8(:,2));
    
    
    % Orthogonalize
    X_GLM_1(:,2) = X_GLM_1(:,2) - (X_GLM_1(:,1)'*X_GLM_1(:,2))*X_GLM_1(:,1);
    X_GLM_2(:,2) = X_GLM_2(:,2) - (X_GLM_2(:,1)'*X_GLM_2(:,2))*X_GLM_2(:,1);
    X_GLM_3(:,2) = X_GLM_3(:,2) - (X_GLM_3(:,1)'*X_GLM_3(:,2))*X_GLM_3(:,1);
    X_GLM_4(:,2) = X_GLM_4(:,2) - (X_GLM_4(:,1)'*X_GLM_4(:,2))*X_GLM_4(:,1);
    X_GLM_5(:,2) = X_GLM_5(:,2) - (X_GLM_5(:,1)'*X_GLM_5(:,2))*X_GLM_5(:,1);
    X_GLM_6(:,2) = X_GLM_6(:,2) - (X_GLM_6(:,1)'*X_GLM_6(:,2))*X_GLM_6(:,1);
    X_GLM_7(:,2) = X_GLM_7(:,2) - (X_GLM_7(:,1)'*X_GLM_7(:,2))*X_GLM_7(:,1);
    X_GLM_8(:,2) = X_GLM_8(:,2) - (X_GLM_8(:,1)'*X_GLM_8(:,2))*X_GLM_8(:,1);
    
    
    xtxxt_GLM_1 = inv(X_GLM_1'*X_GLM_1)*X_GLM_1';
    xtxxt_GLM_2 = inv(X_GLM_2'*X_GLM_2)*X_GLM_2';
    xtxxt_GLM_3 = inv(X_GLM_3'*X_GLM_3)*X_GLM_3';
    xtxxt_GLM_4 = inv(X_GLM_4'*X_GLM_4)*X_GLM_4';
    xtxxt_GLM_5 = inv(X_GLM_5'*X_GLM_5)*X_GLM_5';
    xtxxt_GLM_6 = inv(X_GLM_6'*X_GLM_6)*X_GLM_6';
    xtxxt_GLM_7 = inv(X_GLM_7'*X_GLM_7)*X_GLM_7';
    xtxxt_GLM_8 = inv(X_GLM_8'*X_GLM_8)*X_GLM_8';
    
    c = [1 ; 0];
    
    ctxtx_GLM_1 = c'*inv(X_GLM_1'*X_GLM_1)*c;
    ctxtx_GLM_2 = c'*inv(X_GLM_2'*X_GLM_2)*c;
    ctxtx_GLM_3 = c'*inv(X_GLM_3'*X_GLM_3)*c;
    ctxtx_GLM_4 = c'*inv(X_GLM_4'*X_GLM_4)*c;
    ctxtx_GLM_5 = c'*inv(X_GLM_5'*X_GLM_5)*c;
    ctxtx_GLM_6 = c'*inv(X_GLM_6'*X_GLM_6)*c;
    ctxtx_GLM_7 = c'*inv(X_GLM_7'*X_GLM_7)*c;
    ctxtx_GLM_8 = c'*inv(X_GLM_8'*X_GLM_8)*c;
    
    %----------------
    
    permutation_matrix = zeros(st,number_of_permutations);
    
    for p = 1:number_of_permutations
        permutation_matrix(:,p) = randperm(st);
    end
           
    threshold = icdf('t',1-0.001,st);
    
    
    % Perform permutations on GPU(s)
    
    tic
    [activity_maps_1, activity_maps_2, activity_maps_3, activity_maps_4, activity_maps_5, activity_maps_6, activity_maps_7, activity_maps_8, ...
        voxel_thresholds_1, voxel_thresholds_2, voxel_thresholds_3, voxel_thresholds_4, voxel_thresholds_5, voxel_thresholds_6, voxel_thresholds_7, voxel_thresholds_8 ...
        cluster_thresholds_1, cluster_thresholds_2, cluster_thresholds_3, cluster_thresholds_4, cluster_thresholds_5, cluster_thresholds_6, cluster_thresholds_7, cluster_thresholds_8] ...
        = Preprocessing_and_Permuted_GLM_MultiGPU_Rest_Cluster(fMRI_volumes, brain_voxels, certainty, sum(brain_voxels(:)), f1, f2, f3, number_of_iterations, motion_compensation_filter_size, ...
        smoothing_filters_x, smoothing_filters_y, smoothing_filters_z, alpha_smoothing_filter_x, alpha_smoothing_filter_y, alpha_smoothing_filter_z, smoothing_filter_size, smoothed_certainties, smoothed_alpha_certainty, ...
        X_Detrend, xtxxt_Detrend', ...
        X_GLM_1, xtxxt_GLM_1', ctxtx_GLM_1, X_GLM_2, xtxxt_GLM_2', ctxtx_GLM_2, X_GLM_3, xtxxt_GLM_3', ctxtx_GLM_3, X_GLM_4, xtxxt_GLM_4', ctxtx_GLM_4, ...
        X_GLM_5, xtxxt_GLM_5', ctxtx_GLM_5, X_GLM_6, xtxxt_GLM_6', ctxtx_GLM_6, X_GLM_7, xtxxt_GLM_7', ctxtx_GLM_7, X_GLM_8, xtxxt_GLM_8', ctxtx_GLM_8, ...
        number_of_permutations, uint16((permutation_matrix - 1)), number_of_GPUs, threshold);
    toc
    
        for smoothing = 8:8;
            
            % 1
            % Voxel
            voxel_threshold_1 = voxel_thresholds_1(smoothing+1);
            load voxel_thresholds_rpt_boxcar10.mat
            voxel_thresholds_rpt(file_number,smoothing+1) = voxel_threshold_1;
            save voxel_thresholds_rpt_boxcar10.mat voxel_thresholds_rpt
            
            activity_map_1 = activity_maps_1(:,:,:,smoothing+1);
            significant = sum(activity_map_1(:) > voxel_threshold_1);
            load significant_voxels_rpt_boxcar10.mat
            if significant > 0
                significant_voxels_rpt(file_number,smoothing+1) = 1;
            else
                significant_voxels_rpt(file_number,smoothing+1) = 0;
            end
            save significant_voxels_rpt_boxcar10.mat significant_voxels_rpt
            
            load max_test_values_rpt_boxcar10.mat
            max_test_values_rpt(file_number,smoothing+1) = max(activity_map_1(:));
            save max_test_values_rpt_boxcar10.mat max_test_values_rpt
            
            % Cluster
            cluster_threshold_1 = cluster_thresholds_1(smoothing+1);
            load cluster_thresholds_rpt_boxcar10.mat
            cluster_thresholds_rpt(file_number,smoothing+1) = cluster_threshold_1;
            save cluster_thresholds_rpt_boxcar10.mat cluster_thresholds_rpt
            
            % Get size of largest cluster
            [L,nc] = spm_bwlabel(double(activity_map_1 > threshold),18);
            max_cluster_size = 0;
            for c=1:nc
                summ = sum(L(:) == c);
                if summ > max_cluster_size
                    max_cluster_size = summ;
                end
            end
            
            if max_cluster_size >= cluster_threshold_1
                significant = 1;
            else
                significant = 0;
            end
            
            load significant_clusters_rpt_boxcar10.mat
            significant_clusters_rpt(file_number,smoothing+1) = significant;
            save significant_clusters_rpt_boxcar10.mat significant_clusters_rpt
            
            load max_cluster_sizes_rpt_boxcar10.mat
            max_cluster_sizes_rpt(file_number,smoothing+1) = max_cluster_size;
            save max_cluster_sizes_rpt_boxcar10.mat max_cluster_sizes_rpt
            
            %----------------------------------------------------------------------------------
            
            % 2
            % Voxel
            voxel_threshold_2 = voxel_thresholds_2(smoothing+1);
            load voxel_thresholds_rpt_boxcar15.mat
            voxel_thresholds_rpt(file_number,smoothing+1) = voxel_threshold_2;
            save voxel_thresholds_rpt_boxcar15.mat voxel_thresholds_rpt
            
            activity_map_2 = activity_maps_2(:,:,:,smoothing+1);
            significant = sum(activity_map_2(:) > voxel_threshold_2);
            load significant_voxels_rpt_boxcar15.mat
            if significant > 0
                significant_voxels_rpt(file_number,smoothing+1) = 1;
            else
                significant_voxels_rpt(file_number,smoothing+1) = 0;
            end
            save significant_voxels_rpt_boxcar15.mat significant_voxels_rpt
            
            load max_test_values_rpt_boxcar15.mat
            max_test_values_rpt(file_number,smoothing+1) = max(activity_map_2(:));
            save max_test_values_rpt_boxcar15.mat max_test_values_rpt
            
            % Cluster
            cluster_threshold_2 = cluster_thresholds_2(smoothing+1);
            load cluster_thresholds_rpt_boxcar15.mat
            cluster_thresholds_rpt(file_number,smoothing+1) = cluster_threshold_2;
            save cluster_thresholds_rpt_boxcar15.mat cluster_thresholds_rpt
            
            % Get size of largest cluster
            [L,nc] = spm_bwlabel(double(activity_map_2 > threshold),18);
            max_cluster_size = 0;
            for c=1:nc
                summ = sum(L(:) == c);
                if summ > max_cluster_size
                    max_cluster_size = summ;
                end
            end
            
            if max_cluster_size >= cluster_threshold_2
                significant = 1;
            else
                significant = 0;
            end
            
            load significant_clusters_rpt_boxcar15.mat
            significant_clusters_rpt(file_number,smoothing+1) = significant;
            save significant_clusters_rpt_boxcar15.mat significant_clusters_rpt
            
            load max_cluster_sizes_rpt_boxcar15.mat
            max_cluster_sizes_rpt(file_number,smoothing+1) = max_cluster_size;
            save max_cluster_sizes_rpt_boxcar15.mat max_cluster_sizes_rpt
            
            %----------------------------------------------------------------------------------
            
            % 3
            % Voxel
            voxel_threshold_3 = voxel_thresholds_3(smoothing+1);
            load voxel_thresholds_rpt_boxcar20.mat
            voxel_thresholds_rpt(file_number,smoothing+1) = voxel_threshold_3;
            save voxel_thresholds_rpt_boxcar20.mat voxel_thresholds_rpt
            
            activity_map_3 = activity_maps_3(:,:,:,smoothing+1);
            significant = sum(activity_map_3(:) > voxel_threshold_3);
            load significant_voxels_rpt_boxcar20.mat
            if significant > 0
                significant_voxels_rpt(file_number,smoothing+1) = 1;
            else
                significant_voxels_rpt(file_number,smoothing+1) = 0;
            end
            save significant_voxels_rpt_boxcar20.mat significant_voxels_rpt
            
            load max_test_values_rpt_boxcar20.mat
            max_test_values_rpt(file_number,smoothing+1) = max(activity_map_3(:));
            save max_test_values_rpt_boxcar20.mat max_test_values_rpt
            
            % Cluster
            cluster_threshold_3 = cluster_thresholds_3(smoothing+1);
            load cluster_thresholds_rpt_boxcar20.mat
            cluster_thresholds_rpt(file_number,smoothing+1) = cluster_threshold_3;
            save cluster_thresholds_rpt_boxcar20.mat cluster_thresholds_rpt
            
            % Get size of largest cluster
            [L,nc] = spm_bwlabel(double(activity_map_3 > threshold),18);
            max_cluster_size = 0;
            for c=1:nc
                summ = sum(L(:) == c);
                if summ > max_cluster_size
                    max_cluster_size = summ;
                end
            end
            
            if max_cluster_size >= cluster_threshold_3
                significant = 1;
            else
                significant = 0;
            end
            
            load significant_clusters_rpt_boxcar20.mat
            significant_clusters_rpt(file_number,smoothing+1) = significant;
            save significant_clusters_rpt_boxcar20.mat significant_clusters_rpt
            
            load max_cluster_sizes_rpt_boxcar20.mat
            max_cluster_sizes_rpt(file_number,smoothing+1) = max_cluster_size;
            save max_cluster_sizes_rpt_boxcar20.mat max_cluster_sizes_rpt
            
            %----------------------------------------------------------------------------------
            
            % 4
            % Voxel
            voxel_threshold_4 = voxel_thresholds_4(smoothing+1);
            load voxel_thresholds_rpt_boxcar30.mat
            voxel_thresholds_rpt(file_number,smoothing+1) = voxel_threshold_4;
            save voxel_thresholds_rpt_boxcar30.mat voxel_thresholds_rpt
            
            activity_map_4 = activity_maps_4(:,:,:,smoothing+1);
            significant = sum(activity_map_4(:) > voxel_threshold_4);
            load significant_voxels_rpt_boxcar30.mat
            if significant > 0
                significant_voxels_rpt(file_number,smoothing+1) = 1;
            else
                significant_voxels_rpt(file_number,smoothing+1) = 0;
            end
            save significant_voxels_rpt_boxcar30.mat significant_voxels_rpt
            
            load max_test_values_rpt_boxcar30.mat
            max_test_values_rpt(file_number,smoothing+1) = max(activity_map_4(:));
            save max_test_values_rpt_boxcar30.mat max_test_values_rpt
            
            % Cluster
            cluster_threshold_4 = cluster_thresholds_4(smoothing+1);
            load cluster_thresholds_rpt_boxcar30.mat
            cluster_thresholds_rpt(file_number,smoothing+1) = cluster_threshold_4;
            save cluster_thresholds_rpt_boxcar30.mat cluster_thresholds_rpt
            
            % Get size of largest cluster
            [L,nc] = spm_bwlabel(double(activity_map_4 > threshold),18);
            max_cluster_size = 0;
            for c=1:nc
                summ = sum(L(:) == c);
                if summ > max_cluster_size
                    max_cluster_size = summ;
                end
            end
            
            if max_cluster_size >= cluster_threshold_4
                significant = 1;
            else
                significant = 0;
            end
            
            load significant_clusters_rpt_boxcar30.mat
            significant_clusters_rpt(file_number,smoothing+1) = significant;
            save significant_clusters_rpt_boxcar30.mat significant_clusters_rpt
            
            load max_cluster_sizes_rpt_boxcar30.mat
            max_cluster_sizes_rpt(file_number,smoothing+1) = max_cluster_size;
            save max_cluster_sizes_rpt_boxcar30.mat max_cluster_sizes_rpt
            
            %----------------------------------------------------------------------------------
            
            % 5
            % Voxel
            voxel_threshold_5 = voxel_thresholds_5(smoothing+1);
            load voxel_thresholds_rpt_event1.mat
            voxel_thresholds_rpt(file_number,smoothing+1) = voxel_threshold_5;
            save voxel_thresholds_rpt_event1.mat voxel_thresholds_rpt
            
            activity_map_5 = activity_maps_5(:,:,:,smoothing+1);
            significant = sum(activity_map_5(:) > voxel_threshold_5);
            load significant_voxels_rpt_event1.mat
            if significant > 0
                significant_voxels_rpt(file_number,smoothing+1) = 1;
            else
                significant_voxels_rpt(file_number,smoothing+1) = 0;
            end
            save significant_voxels_rpt_event1.mat significant_voxels_rpt
            
            load max_test_values_rpt_event1.mat
            max_test_values_rpt(file_number,smoothing+1) = max(activity_map_5(:));
            save max_test_values_rpt_event1.mat max_test_values_rpt
            
            % Cluster
            cluster_threshold_5 = cluster_thresholds_5(smoothing+1);
            load cluster_thresholds_rpt_event1.mat
            cluster_thresholds_rpt(file_number,smoothing+1) = cluster_threshold_5;
            save cluster_thresholds_rpt_event1.mat cluster_thresholds_rpt
            
            % Get size of largest cluster
            [L,nc] = spm_bwlabel(double(activity_map_5 > threshold),18);
            max_cluster_size = 0;
            for c=1:nc
                summ = sum(L(:) == c);
                if summ > max_cluster_size
                    max_cluster_size = summ;
                end
            end
            
            if max_cluster_size >= cluster_threshold_5
                significant = 1;
            else
                significant = 0;
            end
            
            load significant_clusters_rpt_event1.mat
            significant_clusters_rpt(file_number,smoothing+1) = significant;
            save significant_clusters_rpt_event1.mat significant_clusters_rpt
            
            load max_cluster_sizes_rpt_event1.mat
            max_cluster_sizes_rpt(file_number,smoothing+1) = max_cluster_size;
            save max_cluster_sizes_rpt_event1.mat max_cluster_sizes_rpt
            
            %----------------------------------------------------------------------------------
            
            % 6
            % Voxel
            voxel_threshold_6 = voxel_thresholds_6(smoothing+1);
            load voxel_thresholds_rpt_event2.mat
            voxel_thresholds_rpt(file_number,smoothing+1) = voxel_threshold_6;
            save voxel_thresholds_rpt_event2.mat voxel_thresholds_rpt
            
            activity_map_6 = activity_maps_6(:,:,:,smoothing+1);
            significant = sum(activity_map_6(:) > voxel_threshold_6);
            load significant_voxels_rpt_event2.mat
            if significant > 0
                significant_voxels_rpt(file_number,smoothing+1) = 1;
            else
                significant_voxels_rpt(file_number,smoothing+1) = 0;
            end
            save significant_voxels_rpt_event2.mat significant_voxels_rpt
            
            load max_test_values_rpt_event2.mat
            max_test_values_rpt(file_number,smoothing+1) = max(activity_map_6(:));
            save max_test_values_rpt_event2.mat max_test_values_rpt
            
            % Cluster
            cluster_threshold_6 = cluster_thresholds_6(smoothing+1);
            load cluster_thresholds_rpt_event2.mat
            cluster_thresholds_rpt(file_number,smoothing+1) = cluster_threshold_6;
            save cluster_thresholds_rpt_event2.mat cluster_thresholds_rpt
            
            % Get size of largest cluster
            [L,nc] = spm_bwlabel(double(activity_map_6 > threshold),18);
            max_cluster_size = 0;
            for c=1:nc
                summ = sum(L(:) == c);
                if summ > max_cluster_size
                    max_cluster_size = summ;
                end
            end
            
            if max_cluster_size >= cluster_threshold_6
                significant = 1;
            else
                significant = 0;
            end
            
            load significant_clusters_rpt_event2.mat
            significant_clusters_rpt(file_number,smoothing+1) = significant;
            save significant_clusters_rpt_event2.mat significant_clusters_rpt
            
            load max_cluster_sizes_rpt_event2.mat
            max_cluster_sizes_rpt(file_number,smoothing+1) = max_cluster_size;
            save max_cluster_sizes_rpt_event2.mat max_cluster_sizes_rpt
            
            %----------------------------------------------------------------------------------
            
             % 7
            % Voxel
            voxel_threshold_7 = voxel_thresholds_7(smoothing+1);
            load voxel_thresholds_rpt_event3.mat
            voxel_thresholds_rpt(file_number,smoothing+1) = voxel_threshold_7;
            save voxel_thresholds_rpt_event3.mat voxel_thresholds_rpt
            
            activity_map_7 = activity_maps_7(:,:,:,smoothing+1);
            significant = sum(activity_map_7(:) > voxel_threshold_7);
            load significant_voxels_rpt_event3.mat
            if significant > 0
                significant_voxels_rpt(file_number,smoothing+1) = 1;
            else
                significant_voxels_rpt(file_number,smoothing+1) = 0;
            end
            save significant_voxels_rpt_event3.mat significant_voxels_rpt
            
            load max_test_values_rpt_event3.mat
            max_test_values_rpt(file_number,smoothing+1) = max(activity_map_7(:));
            save max_test_values_rpt_event3.mat max_test_values_rpt
            
            % Cluster
            cluster_threshold_7 = cluster_thresholds_7(smoothing+1);
            load cluster_thresholds_rpt_event3.mat
            cluster_thresholds_rpt(file_number,smoothing+1) = cluster_threshold_7;
            save cluster_thresholds_rpt_event3.mat cluster_thresholds_rpt
            
            % Get size of largest cluster
            [L,nc] = spm_bwlabel(double(activity_map_7 > threshold),18);
            max_cluster_size = 0;
            for c=1:nc
                summ = sum(L(:) == c);
                if summ > max_cluster_size
                    max_cluster_size = summ;
                end
            end
            
            if max_cluster_size >= cluster_threshold_7
                significant = 1;
            else
                significant = 0;
            end
            
            load significant_clusters_rpt_event3.mat
            significant_clusters_rpt(file_number,smoothing+1) = significant;
            save significant_clusters_rpt_event3.mat significant_clusters_rpt
            
            load max_cluster_sizes_rpt_event3.mat
            max_cluster_sizes_rpt(file_number,smoothing+1) = max_cluster_size;
            save max_cluster_sizes_rpt_event3.mat max_cluster_sizes_rpt
            
            %----------------------------------------------------------------------------------
            
             % 8
            % Voxel
            voxel_threshold_8 = voxel_thresholds_8(smoothing+1);
            load voxel_thresholds_rpt_event4.mat
            voxel_thresholds_rpt(file_number,smoothing+1) = voxel_threshold_8;
            save voxel_thresholds_rpt_event4.mat voxel_thresholds_rpt
            
            activity_map_8 = activity_maps_8(:,:,:,smoothing+1);
            significant = sum(activity_map_8(:) > voxel_threshold_8);
            load significant_voxels_rpt_event4.mat
            if significant > 0
                significant_voxels_rpt(file_number,smoothing+1) = 1;
            else
                significant_voxels_rpt(file_number,smoothing+1) = 0;
            end
            save significant_voxels_rpt_event4.mat significant_voxels_rpt
            
            load max_test_values_rpt_event4.mat
            max_test_values_rpt(file_number,smoothing+1) = max(activity_map_8(:));
            save max_test_values_rpt_event4.mat max_test_values_rpt
            
            % Cluster
            cluster_threshold_8 = cluster_thresholds_8(smoothing+1);
            load cluster_thresholds_rpt_event4.mat
            cluster_thresholds_rpt(file_number,smoothing+1) = cluster_threshold_8;
            save cluster_thresholds_rpt_event4.mat cluster_thresholds_rpt
            
            % Get size of largest cluster
            [L,nc] = spm_bwlabel(double(activity_map_8 > threshold),18);
            max_cluster_size = 0;
            for c=1:nc
                summ = sum(L(:) == c);
                if summ > max_cluster_size
                    max_cluster_size = summ;
                end
            end
            
            if max_cluster_size >= cluster_threshold_8
                significant = 1;
            else
                significant = 0;
            end
            
            load significant_clusters_rpt_event4.mat
            significant_clusters_rpt(file_number,smoothing+1) = significant;
            save significant_clusters_rpt_event4.mat significant_clusters_rpt
            
            load max_cluster_sizes_rpt_event4.mat
            max_cluster_sizes_rpt(file_number,smoothing+1) = max_cluster_size;
            save max_cluster_sizes_rpt_event4.mat max_cluster_sizes_rpt
           
            %----------------------------------------------------------------------------------
            
        end
        
    %end
    
    pause(5)
    
end


