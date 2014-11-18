% voxel_thresholds_spm = zeros(1484,13);
% save voxel_thresholds_spm_boxcar10_files1_375.mat voxel_thresholds_spm
% save voxel_thresholds_spm_boxcar15_files1_375.mat voxel_thresholds_spm
% save voxel_thresholds_spm_boxcar20_files1_375.mat voxel_thresholds_spm
% save voxel_thresholds_spm_boxcar30_files1_375.mat voxel_thresholds_spm
% save voxel_thresholds_spm_event1_files1_375.mat voxel_thresholds_spm
% save voxel_thresholds_spm_event2_files1_375.mat voxel_thresholds_spm
% save voxel_thresholds_spm_event3_files1_375.mat voxel_thresholds_spm
% save voxel_thresholds_spm_event4_files1_375.mat voxel_thresholds_spm
% 
% 
% cluster_thresholds_spm = zeros(1484,13);
% save cluster_thresholds_spm_boxcar10_files1_375.mat cluster_thresholds_spm
% save cluster_thresholds_spm_boxcar15_files1_375.mat cluster_thresholds_spm
% save cluster_thresholds_spm_boxcar20_files1_375.mat cluster_thresholds_spm
% save cluster_thresholds_spm_boxcar30_files1_375.mat cluster_thresholds_spm
% save cluster_thresholds_spm_event1_files1_375.mat cluster_thresholds_spm
% save cluster_thresholds_spm_event2_files1_375.mat cluster_thresholds_spm
% save cluster_thresholds_spm_event3_files1_375.mat cluster_thresholds_spm
% save cluster_thresholds_spm_event4_files1_375.mat cluster_thresholds_spm
% 
% max_test_values_spm = zeros(1484,13);
% save max_test_values_spm_boxcar10_files1_375.mat max_test_values_spm
% save max_test_values_spm_boxcar15_files1_375.mat max_test_values_spm
% save max_test_values_spm_boxcar20_files1_375.mat max_test_values_spm
% save max_test_values_spm_boxcar30_files1_375.mat max_test_values_spm
% save max_test_values_spm_event1_files1_375.mat max_test_values_spm
% save max_test_values_spm_event2_files1_375.mat max_test_values_spm
% save max_test_values_spm_event3_files1_375.mat max_test_values_spm
% save max_test_values_spm_event4_files1_375.mat max_test_values_spm
% 
% max_cluster_sizes_spm = zeros(1484,13);
% save max_cluster_sizes_spm_boxcar10_files1_375.mat max_cluster_sizes_spm
% save max_cluster_sizes_spm_boxcar15_files1_375.mat max_cluster_sizes_spm
% save max_cluster_sizes_spm_boxcar20_files1_375.mat max_cluster_sizes_spm
% save max_cluster_sizes_spm_boxcar30_files1_375.mat max_cluster_sizes_spm
% save max_cluster_sizes_spm_event1_files1_375.mat max_cluster_sizes_spm
% save max_cluster_sizes_spm_event2_files1_375.mat max_cluster_sizes_spm
% save max_cluster_sizes_spm_event3_files1_375.mat max_cluster_sizes_spm
% save max_cluster_sizes_spm_event4_files1_375.mat max_cluster_sizes_spm
% 
% significant_voxels_spm = zeros(1484,13);
% save significant_voxels_spm_boxcar10_files1_375.mat significant_voxels_spm
% save significant_voxels_spm_boxcar15_files1_375.mat significant_voxels_spm
% save significant_voxels_spm_boxcar20_files1_375.mat significant_voxels_spm
% save significant_voxels_spm_boxcar30_files1_375.mat significant_voxels_spm
% save significant_voxels_spm_event1_files1_375.mat significant_voxels_spm
% save significant_voxels_spm_event2_files1_375.mat significant_voxels_spm
% save significant_voxels_spm_event3_files1_375.mat significant_voxels_spm
% save significant_voxels_spm_event4_files1_375.mat significant_voxels_spm
% 
% significant_clusters_spm = zeros(1484,13);
% save significant_clusters_spm_boxcar10_files1_375.mat significant_clusters_spm
% save significant_clusters_spm_boxcar15_files1_375.mat significant_clusters_spm
% save significant_clusters_spm_boxcar20_files1_375.mat significant_clusters_spm
% save significant_clusters_spm_boxcar30_files1_375.mat significant_clusters_spm
% save significant_clusters_spm_event1_files1_375.mat significant_clusters_spm
% save significant_clusters_spm_event2_files1_375.mat significant_clusters_spm
% save significant_clusters_spm_event3_files1_375.mat significant_clusters_spm
% save significant_clusters_spm_event4_files1_375.mat significant_clusters_spm

%%

clear all
clc
close all

load rand_events1.mat
onsets1 = onsets;
durations1 = durations;

load rand_events2.mat
onsets2 = onsets;
durations2 = durations;

addpath('E:\rest_fMRI');
addpath('E:\spm8');

for file_number = 1:375
    
    tic
    
    if file_number < 10
        filename = ['E:\rest_fMRI\func\subject_000' num2str(file_number) '\rest_' num2str(file_number) '.nii'];
        filedirectory = ['E:\rest_fMRI\func\subject_000' num2str(file_number) '\'];
    elseif file_number < 100
        filename = ['E:\rest_fMRI\func\subject_00' num2str(file_number) '\rest_' num2str(file_number) '.nii'];
        filedirectory = ['E:\rest_fMRI\func\subject_00' num2str(file_number) '\'];
    elseif file_number < 1000
        filename = ['E:\rest_fMRI\func\subject_0' num2str(file_number) '\rest_' num2str(file_number) '.nii'];
        filedirectory = ['E:\rest_fMRI\func\subject_0' num2str(file_number) '\'];
    else
        filename = ['E:\rest_fMRI\func\subject_' num2str(file_number) '\rest_' num2str(file_number) '.nii'];
        filedirectory = ['E:\rest_fMRI\func\subject_' num2str(file_number) '\'];
    end
    
    V = spm_vol(filename);
    
    % Get number of time points
    st = size(V,1);
    st
    
    % Get repetition time
    a = V.private;
    TR = a.timing.tspace;
    
    
    %% Path containing data
    %--------------------------------------------------------------------------
    data_path = 'E:\rest_fMRI\func_only_1\';
    data_path2 = 'E:\rest_fMRI\';
    
    %% Initialise SPM defaults
    %--------------------------------------------------------------------------
    spm('Defaults','fMRI');
    
    spm_jobman('initcfg'); % useful in SPM8 only
    
    
    %% WORKING DIRECTORY (useful for .ps only)
    %--------------------------------------------------------------------------
    clear pjobs
    pjobs{1}.util{1}.cdir.directory = cellstr(data_path);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPATIAL PREPROCESSING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Select functional and structural scans
    %--------------------------------------------------------------------------   
    
    % Select data
    if file_number < 10
        f = spm_select('FPList', fullfile(data_path2,['func\subject_000' num2str(file_number) '\']),'func*') ;
    elseif file_number < 100
        f = spm_select('FPList', fullfile(data_path2,['func\subject_00' num2str(file_number) '\']),'func*') ;
    elseif file_number < 1000
        f = spm_select('FPList', fullfile(data_path2,['func\subject_0' num2str(file_number) '\']),'func*') ;
    else
        f = spm_select('FPList', fullfile(data_path2,['func\subject_' num2str(file_number) '\']),'func*') ;
    end
    % Comment if preprocessed files needs to be calculated
    f = f(1:st,:);
    
    
    %% REALIGN
    %--------------------------------------------------------------------------
    pjobs{2}.spatial{1}.realign{1}.estwrite.data{1} = cellstr(f);
       
    %% SMOOTHING
    %--------------------------------------------------------------------------
    
    %smoothings = 4:16;
    smoothings = 4:2:16;
    %for smoothing = 1:13
    for smoothing = 1:7    
        pjobs{2}.spatial{2 + smoothing - 1}.smooth.data = editfilenames(f,'prefix','r');
        pjobs{2}.spatial{2 + smoothing - 1}.smooth.fwhm = [smoothings(smoothing) smoothings(smoothing) smoothings(smoothing)];
        pjobs{2}.spatial{2 + smoothing - 1}.smooth.prefix = ['s' num2str(smoothings(smoothing))];
    end
    smoothing_saves = [1 3 5 7 9 11 13];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FINISH
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save('rest_batch_preprocessing.mat','pjobs');
    %spm_jobman('interactive',jobs); % open a GUI containing all the setup
    error1 = 0;
    try
        %spm_jobman('run',pjobs);        % run preprocessing
    catch err
        err
        error1 = 1;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CLASSICAL STATISTICAL ANALYSIS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if error1 == 0
        
        for exp = 1:8
            %for smoothing = 1:13
            for smoothing = 1:7
            smoothing_save = smoothing_saves(smoothing);
                
                clear jobs
                jobs{1}.util{1}.cdir.directory = cellstr(data_path);
                
                % Remove old SPM.mat
                spm_file = fullfile('E:\rest_fMRI\func_only_1\classical','SPM.mat');
                if exist(spm_file,'file')==2
                    %system(['rm' spm_file]); % Linux
                    delete('E:\rest_fMRI\func_only_1\classical\SPM.mat')
                end
                
                %% OUTPUT DIRECTORY
                %--------------------------------------------------------------------------
                jobs{1}.util{1}.md.basedir = cellstr(data_path);
                jobs{1}.util{1}.md.name = 'classical';
                
                %% MODEL SPECIFICATION AND ESTIMATION
                %--------------------------------------------------------------------------
                jobs{2}.stats{1}.fmri_spec.dir = cellstr(fullfile(data_path,'classical'));
                jobs{2}.stats{1}.fmri_spec.timing.units = 'secs';
                jobs{2}.stats{1}.fmri_spec.timing.RT = TR;
                jobs{2}.stats{1}.fmri_spec.sess.scans = editfilenames(f,'prefix',['s' num2str(smoothings(smoothing)) 'r']);
                jobs{2}.stats{1}.fmri_spec.sess.cond.name = 'active';
%                 
                if exp == 1  % 10 s boxcar
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = 10:20:(st*TR);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = 10;
                elseif exp == 2 % 15 s boxcar
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = 15:30:(st*TR);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = 15;
                elseif exp == 3 % 20 s boxcar
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = 20:40:(st*TR);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = 20;
                elseif exp == 4 % 30 s boxcar
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = 30:60:(st*TR);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = 30;
                elseif exp == 5 % event 1 (2, 6)
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = 8:8:(st*TR);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = 2;
                elseif exp == 6 % event 2 (4, 8)
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = 12:12:(st*TR);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = 4;
                elseif exp == 7 % event 3 (rand 1)
                    
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
                    
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = onsets1(1:last);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = durations1(1:last);

                elseif exp == 8 % event 4 (rand 2)
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
                    
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = onsets2(1:last);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = durations2(1:last);
                end
                
                % Use temporal derivatives
                jobs{2}.stats{1}.fmri_spec.bases.hrf.derivs = [1 0];
                
                % Use global normalization or not
                jobs{2}.stats{1}.fmri_spec.global = 'Scaling';
                %jobs{2}.stats{1}.fmri_spec.global = 'None';
                
                % Use motion regressors in design matrix or not
                %jobs{2}.stats{1}.fmri_spec.sess.multi_reg = {[filedirectory 'rp_func_001.txt']};
                
                jobs{2}.stats{2}.fmri_est.spmmat = cellstr(fullfile(data_path,'classical','SPM.mat'));
                
                %% INFERENCE
                %--------------------------------------------------------------------------
                jobs{2}.stats{3}.con.spmmat = cellstr(fullfile(data_path,'classical','SPM.mat'));
                jobs{2}.stats{3}.con.consess{1}.tcon = struct('name','active > rest','convec', 1,'sessrep','none');
                
                jobs{2}.stats{4}.results.spmmat = cellstr(fullfile(data_path,'classical','SPM.mat'));
                jobs{2}.stats{4}.results.conspec.contrasts = Inf;
                jobs{2}.stats{4}.results.conspec.threshdesc = 'none';
                jobs{2}.stats{4}.results.conspec.thresh = 0.001; % for cluster based threshold
                jobs{2}.stats{4}.results.conspec.extent = 0;                                        
                                                              
                save('rest_batch_analysis.mat','jobs');
                
                error2 = 0;
                try
                    spm_jobman('run',jobs);        % run analysis
                catch err
                    err
                    error2 = 1;
                end
                
                if error2 == 0
                    
                    if exp == 1                                                
                        
                        % Get corrected threshold
                        corrected_threshold = xSPM.uc(1)
                        % Save corrected threshold
                        load voxel_thresholds_spm_boxcar10_files1_375.mat
                        voxel_thresholds_spm(file_number,smoothing_save) = corrected_threshold;
                        save voxel_thresholds_spm_boxcar10_files1_375.mat voxel_thresholds_spm
                        
                        % Get t-map
                        V = spm_vol('spmT_0001.hdr');
                        [tmap,aa] = spm_read_vols(V);
                        
                        % Check the number of significantly active voxels
                        significant = sum(tmap(:) > corrected_threshold)
                        
                        % Save if there were significant voxels or not
                        load significant_voxels_spm_boxcar10_files1_375.mat
                        if significant > 0
                            significant_voxels_spm(file_number,smoothing_save) = 1;
                        else
                            significant_voxels_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_voxels_spm_boxcar10_files1_375.mat significant_voxels_spm
                        
                        load max_test_values_spm_boxcar10_files1_375.mat
                        max(tmap(:))
                        max_test_values_spm(file_number,smoothing_save) = max(tmap(:));
                        save max_test_values_spm_boxcar10_files1_375.mat max_test_values_spm
                                
                        [k,Pc] = corrclusth(SPM,0.001,0.05,1:100000);
                        k
                        load cluster_thresholds_spm_boxcar10_files1_375.mat
                        cluster_thresholds_spm(file_number,smoothing_save) = k;
                        save cluster_thresholds_spm_boxcar10_files1_375.mat cluster_thresholds_spm
               
                        indices = find(tmap>xSPM.u);
                        max_cluster = max_extent(tmap, indices)                       
                        load max_cluster_sizes_spm_boxcar10_files1_375.mat
                        max_cluster_sizes_spm(file_number,smoothing_save) = max_cluster;
                        save max_cluster_sizes_spm_boxcar10_files1_375.mat max_cluster_sizes_spm
                                                                       
                        load significant_clusters_spm_boxcar10_files1_375.mat
                        if max_cluster >= k
                            significant_clusters_spm(file_number,smoothing_save) = 1;
                        else
                            significant_clusters_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_clusters_spm_boxcar10_files1_375.mat significant_clusters_spm
                        
                        
                    elseif exp == 2
                        
                        % Get corrected threshold
                        corrected_threshold = xSPM.uc(1)
                        % Save corrected threshold
                        load voxel_thresholds_spm_boxcar15_files1_375.mat
                        voxel_thresholds_spm(file_number,smoothing_save) = corrected_threshold;
                        save voxel_thresholds_spm_boxcar15_files1_375.mat voxel_thresholds_spm
                        
                        % Get t-map
                        V = spm_vol('spmT_0001.hdr');
                        [tmap,aa] = spm_read_vols(V);
                        
                        % Check the number of significantly active voxels
                        significant = sum(tmap(:) > corrected_threshold)
                        
                        % Save if there were significant voxels or not
                        load significant_voxels_spm_boxcar15_files1_375.mat
                        if significant > 0
                            significant_voxels_spm(file_number,smoothing_save) = 1;
                        else
                            significant_voxels_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_voxels_spm_boxcar15_files1_375.mat significant_voxels_spm
                        
                        load max_test_values_spm_boxcar15_files1_375.mat
                        max(tmap(:))
                        max_test_values_spm(file_number,smoothing_save) = max(tmap(:));
                        save max_test_values_spm_boxcar15_files1_375.mat max_test_values_spm
                        
                        [k,Pc] = corrclusth(SPM,0.001,0.05,1:100000);
                        k
                        load cluster_thresholds_spm_boxcar15_files1_375.mat
                        cluster_thresholds_spm(file_number,smoothing_save) = k;
                        save cluster_thresholds_spm_boxcar15_files1_375.mat cluster_thresholds_spm
               
                        indices = find(tmap>xSPM.u);
                        max_cluster = max_extent(tmap, indices)                        
                        load max_cluster_sizes_spm_boxcar15_files1_375.mat
                        max_cluster_sizes_spm(file_number,smoothing_save) = max_cluster;
                        save max_cluster_sizes_spm_boxcar15_files1_375.mat max_cluster_sizes_spm                                                
                        
                        load significant_clusters_spm_boxcar15_files1_375.mat
                        if max_cluster >= k
                            significant_clusters_spm(file_number,smoothing_save) = 1;
                        else
                            significant_clusters_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_clusters_spm_boxcar15_files1_375.mat significant_clusters_spm
                        
                        
                    elseif exp == 3
                        
                        % Get corrected threshold
                        corrected_threshold = xSPM.uc(1)
                        % Save corrected threshold
                        load voxel_thresholds_spm_boxcar20_files1_375.mat
                        voxel_thresholds_spm(file_number,smoothing_save) = corrected_threshold;
                        save voxel_thresholds_spm_boxcar20_files1_375.mat voxel_thresholds_spm
                         
                        % Get t-map
                        V = spm_vol('spmT_0001.hdr');
                        [tmap,aa] = spm_read_vols(V);
                        
                        % Check the number of significantly active voxels
                        significant = sum(tmap(:) > corrected_threshold)
                        
                        % Save if there were significant voxels or not
                        load significant_voxels_spm_boxcar20_files1_375.mat
                        if significant > 0
                            significant_voxels_spm(file_number,smoothing_save) = 1;
                        else
                            significant_voxels_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_voxels_spm_boxcar20_files1_375.mat significant_voxels_spm
                        
                        load max_test_values_spm_boxcar20_files1_375.mat
                        max(tmap(:))
                        max_test_values_spm(file_number,smoothing_save) = max(tmap(:));
                        save max_test_values_spm_boxcar20_files1_375.mat max_test_values_spm
                        
                        [k,Pc] = corrclusth(SPM,0.001,0.05,1:100000);
                        k
                        load cluster_thresholds_spm_boxcar20_files1_375.mat
                        cluster_thresholds_spm(file_number,smoothing_save) = k;
                        save cluster_thresholds_spm_boxcar20_files1_375.mat cluster_thresholds_spm
               
                        indices = find(tmap>xSPM.u);
                        max_cluster = max_extent(tmap, indices)                        
                        load max_cluster_sizes_spm_boxcar20_files1_375.mat
                        max_cluster_sizes_spm(file_number,smoothing_save) = max_cluster;
                        save max_cluster_sizes_spm_boxcar20_files1_375.mat max_cluster_sizes_spm
                                                
                        load significant_clusters_spm_boxcar20_files1_375.mat
                        if max_cluster >= k
                            significant_clusters_spm(file_number,smoothing_save) = 1;
                        else
                            significant_clusters_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_clusters_spm_boxcar20_files1_375.mat significant_clusters_spm
                        
                        
                    elseif exp == 4
                        
                        % Get corrected threshold
                        corrected_threshold = xSPM.uc(1)
                        % Save corrected threshold
                        load voxel_thresholds_spm_boxcar30_files1_375.mat
                        voxel_thresholds_spm(file_number,smoothing_save) = corrected_threshold;
                        save voxel_thresholds_spm_boxcar30_files1_375.mat voxel_thresholds_spm
                        
                        % Get t-map
                        V = spm_vol('spmT_0001.hdr');
                        [tmap,aa] = spm_read_vols(V);
                        
                        % Check the number of significantly active voxels
                        significant = sum(tmap(:) > corrected_threshold)
                        
                        % Save if there were significant voxels or not
                        load significant_voxels_spm_boxcar30_files1_375.mat
                        if significant > 0
                            significant_voxels_spm(file_number,smoothing_save) = 1;
                        else
                            significant_voxels_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_voxels_spm_boxcar30_files1_375.mat significant_voxels_spm
                        
                        load max_test_values_spm_boxcar30_files1_375.mat
                        max(tmap(:))
                        max_test_values_spm(file_number,smoothing_save) = max(tmap(:));
                        save max_test_values_spm_boxcar30_files1_375.mat max_test_values_spm
                        
                        [k,Pc] = corrclusth(SPM,0.001,0.05,1:100000);
                        k
                        load cluster_thresholds_spm_boxcar30_files1_375.mat
                        cluster_thresholds_spm(file_number,smoothing_save) = k;
                        save cluster_thresholds_spm_boxcar30_files1_375.mat cluster_thresholds_spm
               
                        
                        indices = find(tmap>xSPM.u);
                        max_cluster = max_extent(tmap, indices)                        
                        load max_cluster_sizes_spm_boxcar30_files1_375.mat
                        max_cluster_sizes_spm(file_number,smoothing_save) = max_cluster;
                        save max_cluster_sizes_spm_boxcar30_files1_375.mat max_cluster_sizes_spm
                                                                      
                        load significant_clusters_spm_boxcar30_files1_375.mat
                        if max_cluster >= k
                            significant_clusters_spm(file_number,smoothing_save) = 1;
                        else
                            significant_clusters_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_clusters_spm_boxcar30_files1_375.mat significant_clusters_spm
                        
                        
                    elseif exp == 5
                        
                        % Get corrected threshold
                        corrected_threshold = xSPM.uc(1)
                        % Save corrected threshold
                        load voxel_thresholds_spm_event1_files1_375.mat
                        voxel_thresholds_spm(file_number,smoothing_save) = corrected_threshold;
                        save voxel_thresholds_spm_event1_files1_375.mat voxel_thresholds_spm
                         
                        % Get t-map
                        V = spm_vol('spmT_0001.hdr');
                        [tmap,aa] = spm_read_vols(V);
                        
                        % Check the number of significantly active voxels
                        significant = sum(tmap(:) > corrected_threshold)
                        
                        % Save if there were significant voxels or not
                        load significant_voxels_spm_event1_files1_375.mat
                        if significant > 0
                            significant_voxels_spm(file_number,smoothing_save) = 1;
                        else
                            significant_voxels_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_voxels_spm_event1_files1_375.mat significant_voxels_spm
                        
                        load max_test_values_spm_event1_files1_375.mat
                        max(tmap(:))
                        max_test_values_spm(file_number,smoothing_save) = max(tmap(:));
                        save max_test_values_spm_event1_files1_375.mat max_test_values_spm
                        
                        [k,Pc] = corrclusth(SPM,0.001,0.05,1:100000);
                        k
                        load cluster_thresholds_spm_event1_files1_375.mat
                        cluster_thresholds_spm(file_number,smoothing_save) = k;
                        save cluster_thresholds_spm_event1_files1_375.mat cluster_thresholds_spm
               
                        indices = find(tmap>xSPM.u);
                        max_cluster = max_extent(tmap, indices)                        
                        load max_cluster_sizes_spm_event1_files1_375.mat
                        max_cluster_sizes_spm(file_number,smoothing_save) = max_cluster;
                        save max_cluster_sizes_spm_event1_files1_375.mat max_cluster_sizes_spm
                                                
                        load significant_clusters_spm_event1_files1_375.mat
                        if max_cluster >= k
                            significant_clusters_spm(file_number,smoothing_save) = 1;
                        else
                            significant_clusters_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_clusters_spm_event1_files1_375.mat significant_clusters_spm
                                                                        
                    elseif exp == 6
                        
                        % Get corrected threshold
                        corrected_threshold = xSPM.uc(1)
                        load voxel_thresholds_spm_event2_files1_375.mat
                        voxel_thresholds_spm(file_number,smoothing_save) = corrected_threshold;
                        save voxel_thresholds_spm_event2_files1_375.mat voxel_thresholds_spm
                        
                        % Get t-map
                        V = spm_vol('spmT_0001.hdr');
                        [tmap,aa] = spm_read_vols(V);
                        
                        % Check the number of significantly active voxels
                        significant = sum(tmap(:) > corrected_threshold)
                        
                        % Save if there were significant voxels or not
                        load significant_voxels_spm_event2_files1_375.mat
                        if significant > 0
                            significant_voxels_spm(file_number,smoothing_save) = 1;
                        else
                            significant_voxels_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_voxels_spm_event2_files1_375.mat significant_voxels_spm
                        
                        load max_test_values_spm_event2_files1_375.mat
                        max(tmap(:))
                        max_test_values_spm(file_number,smoothing_save) = max(tmap(:));
                        save max_test_values_spm_event2_files1_375.mat max_test_values_spm
                        
                        [k,Pc] = corrclusth(SPM,0.001,0.05,1:100000);
                        k
                        load cluster_thresholds_spm_event2_files1_375.mat
                        cluster_thresholds_spm(file_number,smoothing_save) = k;
                        save cluster_thresholds_spm_event2_files1_375.mat cluster_thresholds_spm
                        
                        indices = find(tmap>xSPM.u);
                        max_cluster = max_extent(tmap, indices)                        
                        load max_cluster_sizes_spm_event2_files1_375.mat
                        max_cluster_sizes_spm(file_number,smoothing_save) = max_cluster;
                        save max_cluster_sizes_spm_event2_files1_375.mat max_cluster_sizes_spm
                                                
                        load significant_clusters_spm_event2_files1_375.mat
                        if max_cluster >= k
                            significant_clusters_spm(file_number,smoothing_save) = 1;
                        else
                            significant_clusters_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_clusters_spm_event2_files1_375.mat significant_clusters_spm
                        
                    elseif exp == 7
                        
                        % Get corrected threshold
                        corrected_threshold = xSPM.uc(1)
                        load voxel_thresholds_spm_event3_files1_375.mat
                        voxel_thresholds_spm(file_number,smoothing_save) = corrected_threshold;
                        save voxel_thresholds_spm_event3_files1_375.mat voxel_thresholds_spm
                        
                        % Get t-map
                        V = spm_vol('spmT_0001.hdr');
                        [tmap,aa] = spm_read_vols(V);
                        
                        % Check the number of significantly active voxels
                        significant = sum(tmap(:) > corrected_threshold)
                        
                        % Save if there were significant voxels or not
                        load significant_voxels_spm_event3_files1_375.mat
                        if significant > 0
                            significant_voxels_spm(file_number,smoothing_save) = 1;
                        else
                            significant_voxels_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_voxels_spm_event3_files1_375.mat significant_voxels_spm
                        
                        load max_test_values_spm_event3_files1_375.mat
                        max(tmap(:))
                        max_test_values_spm(file_number,smoothing_save) = max(tmap(:));
                        save max_test_values_spm_event3_files1_375.mat max_test_values_spm
                        
                        [k,Pc] = corrclusth(SPM,0.001,0.05,1:100000);
                        k
                        load cluster_thresholds_spm_event3_files1_375.mat
                        cluster_thresholds_spm(file_number,smoothing_save) = k;
                        save cluster_thresholds_spm_event3_files1_375.mat cluster_thresholds_spm
                        
                        indices = find(tmap>xSPM.u);
                        max_cluster = max_extent(tmap, indices)                        
                        load max_cluster_sizes_spm_event3_files1_375.mat
                        max_cluster_sizes_spm(file_number,smoothing_save) = max_cluster;
                        save max_cluster_sizes_spm_event3_files1_375.mat max_cluster_sizes_spm
                                                
                        load significant_clusters_spm_event3_files1_375.mat
                        if max_cluster >= k
                            significant_clusters_spm(file_number,smoothing_save) = 1;
                        else
                            significant_clusters_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_clusters_spm_event3_files1_375.mat significant_clusters_spm
                        
                        
                    elseif exp == 8
                        
                        % Get corrected threshold
                        corrected_threshold = xSPM.uc(1)
                        load voxel_thresholds_spm_event4_files1_375.mat
                        voxel_thresholds_spm(file_number,smoothing_save) = corrected_threshold;
                        save voxel_thresholds_spm_event4_files1_375.mat voxel_thresholds_spm
                        
                        % Get t-map
                        V = spm_vol('spmT_0001.hdr');
                        [tmap,aa] = spm_read_vols(V);
                        
                        % Check the number of significantly active voxels
                        significant = sum(tmap(:) > corrected_threshold)
                        
                        % Save if there were significant voxels or not
                        load significant_voxels_spm_event4_files1_375.mat
                        if significant > 0
                            significant_voxels_spm(file_number,smoothing_save) = 1;
                        else
                            significant_voxels_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_voxels_spm_event4_files1_375.mat significant_voxels_spm
                        
                        load max_test_values_spm_event4_files1_375.mat
                        max(tmap(:))
                        max_test_values_spm(file_number,smoothing_save) = max(tmap(:));
                        save max_test_values_spm_event4_files1_375.mat max_test_values_spm
                        
                        indices = find(tmap>xSPM.u);
                        max_cluster = max_extent(tmap, indices)                        
                        load max_cluster_sizes_spm_event4_files1_375.mat
                        max_cluster_sizes_spm(file_number,smoothing_save) = max_cluster;
                        save max_cluster_sizes_spm_event4_files1_375.mat max_cluster_sizes_spm
                        
                        [k,Pc] = corrclusth(SPM,0.001,0.05,1:100000);
                        k
                        load cluster_thresholds_spm_event4_files1_375.mat
                        cluster_thresholds_spm(file_number,smoothing_save) = k;
                        save cluster_thresholds_spm_event4_files1_375.mat cluster_thresholds_spm
                                                
                        load significant_clusters_spm_event4_files1_375.mat
                        if max_cluster >= k
                            significant_clusters_spm(file_number,smoothing_save) = 1;
                        else
                            significant_clusters_spm(file_number,smoothing_save) = 0;
                        end
                        save significant_clusters_spm_event4_files1_375.mat significant_clusters_spm
                        
                    end
                    
                    
                    
            
                    
                end
            end
            
        end
        
    end
    
    file_number
    
    
    toc
    
end








