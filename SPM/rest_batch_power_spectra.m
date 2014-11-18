% 
% 
% TIMESERIES = zeros(512,3);
% save power_spectra_boxcar10.mat TIMESERIES
% save power_spectra_boxcar15.mat TIMESERIES
% save power_spectra_boxcar20.mat TIMESERIES
% save power_spectra_boxcar30.mat TIMESERIES
% 
% save power_spectra_event1.mat TIMESERIES
% save power_spectra_event2.mat TIMESERIES
% save power_spectra_event3.mat TIMESERIES
% save power_spectra_event4.mat TIMESERIES
% 
% rest_TRs = zeros(1484,1);
% save rest_TRs.mat rest_TRs;
% 
% number_of_datasets_TR = zeros(3,8);
% save number_of_datasets_TR.mat number_of_datasets_TR
% 
% mean_residual_variances = zeros(1484,8);
% save mean_residual_variances.mat mean_residual_variances
% 
% mean_spectras = zeros(1484,8,512);
% save mean_spectras.mat mean_spectras


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
addpath('E:\spm8_residuals');

for file_number = 1:1484
    
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
    TR_round = round(TR);
    %load rest_TRs.mat
    %rest_TRs(file_number) = TR;
    %save rest_TRs.mat rest_TRs
    
    % Only use datasets with approximate TR of 1, 2 or 3 seconds
    if (abs(TR - TR_round) < 0.05)
              
        %% Path containing data
        %--------------------------------------------------------------------------
        data_path = 'E:\rest_fMRI\func_only_3\';
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
        f = f(1:st,:);
        
        %% REALIGN
        %--------------------------------------------------------------------------
        pjobs{2}.spatial{1}.realign{1}.estwrite.data{1} = cellstr(f);
        
        %% SMOOTHING
        %--------------------------------------------------------------------------
        
        smoothing = 8;
        pjobs{2}.spatial{2}.smooth.data = editfilenames(f,'prefix','r');
        pjobs{2}.spatial{2}.smooth.fwhm = [smoothing smoothing smoothing];
        pjobs{2}.spatial{2}.smooth.prefix = ['s' num2str(smoothing)];
        
        
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
            
            for exp = 1:1
                
                clear jobs
                jobs{1}.util{1}.cdir.directory = cellstr(data_path);
                
                % Remove old SPM.mat
                spm_file = fullfile('E:\rest_fMRI\func_only_3\classical','SPM.mat');
                if exist(spm_file,'file')==2
                    delete('E:\rest_fMRI\func_only_3\classical\SPM.mat')
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
                jobs{2}.stats{1}.fmri_spec.sess.scans = editfilenames(f,'prefix',['s' num2str(smoothing) 'r']);
                jobs{2}.stats{1}.fmri_spec.sess.cond.name = 'active';
                
                if exp == 1  % 10 s boxcar
                    jobs{2}.stats{1}.fmri_spec.sess.cond.onset = 10:20:(st*TR);
                    jobs{2}.stats{1}.fmri_spec.sess.cond.duration = 10;
                end
                
                % Use temporal derivatives
                jobs{2}.stats{1}.fmri_spec.bases.hrf.derivs = [1 0];
                
                % Use global normalization or not
                jobs{2}.stats{1}.fmri_spec.global = 'Scaling';
                %jobs{2}.stats{1}.fmri_spec.global = 'None';
                
                % Use motion regressors or not
                jobs{2}.stats{1}.fmri_spec.sess.multi_reg = {[filedirectory 'rp_func_001.txt']};
                
                % Use whitening or not
                %jobs{2}.stats{1}.fmri_spec.cvi = 'none';
                jobs{2}.stats{1}.fmri_spec.cvi = 'AR(1)';
                
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
                                                
                        V = spm_vol('spmT_0001.hdr');
                        [tmap,aa] = spm_read_vols(V);
                        % Get size of volume
                        [sy sx sz] = size(tmap);
                        % Get brain mask
                        brain_mask = double(tmap ~= 0);
                        
                        residual_volumes = zeros(sy,sx,sz,st);
                        
                        % Read residual volumes                        
                        for t = 1:st
                            if t < 10
                                V = spm_vol(['ResI_000' num2str(t) '.hdr']);
                            elseif t < 100
                                V = spm_vol(['ResI_00' num2str(t) '.hdr']);
                            elseif t < 1000
                                V = spm_vol(['ResI_0' num2str(t) '.hdr']);
                            end                            
                            [residual_volumes(:,:,:,t),aa] = spm_read_vols(V);
                        end
                        
                        
                        % Calculate Fourier transform of each timeseries
                        % in the mask                      
                        TIMESERIES_ = zeros(512,1);
                        for x = 1:sx
                            for y = 1:sy
                                for z = 1:sz
                                    if (brain_mask(y,x,z) == 1)
                                        % 512 point FFT
                                        timeseries = squeeze(residual_volumes(y,x,z,:));             % Variance 1, standardized                           
                                        timeseries = fft(timeseries,512);
                                        % Calculate power of spectra
                                        TIMESERIES_ = TIMESERIES_ + ((abs(timeseries)).^2)/st;
                                    end
                                end
                            end
                        end
                        % Calculate mean spectra
                        TIMESERIES_ = TIMESERIES_ / sum(brain_mask(:));
                        
                        load mean_spectras.mat
                        mean_spectras(file_number,exp,:) = TIMESERIES_;
                        save mean_spectras.mat mean_spectras
                                                
                        % Remove residual files
                        j = spm_select('List',SPM.swd,'^ResI_.{4}\..{3}$');
                        for  k = 1:size(j,1)
                            spm_unlink(deblank(j(k,:)));
                        end
                                            
                    end
                                                                                                    
                end
            end
            
        end
        
    else
        file_number
        disp('Skipped dataset since TR is not integer')  
    end
    
    file_number
    
    
    toc
    
end








