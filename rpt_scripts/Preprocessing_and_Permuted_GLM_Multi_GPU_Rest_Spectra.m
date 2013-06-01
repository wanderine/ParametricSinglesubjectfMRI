
clear all
close all
clc

%mean_original_spectras = zeros(1484,8,512);
%save mean_original_spectras.mat mean_original_spectras
% 
%mean_whitened_spectras = zeros(1484,8,512);
%save mean_whitened_spectras.mat mean_whitened_spectras


addpath('/home/andek/Testsaker/spm8');
addpath('/home/andek/Testsaker/nifti_matlab');
number_of_permutations = 99;
alpha_FWHM = 8.01;
T_FWHM = 8.01;
smoothing_filter_size = 9;
plott = 1;

smoothing_filters_x = zeros(smoothing_filter_size,16);
smoothing_filters_y = zeros(smoothing_filter_size,16);
smoothing_filters_z = zeros(smoothing_filter_size,16);

for file_number = 1:1484
    file_number
    
    % Load data
    filename = ['/flush/andek/rest_fMRI/func_nii/rest_' num2str(file_number) '.nii'];
    V = spm_vol(filename);
    
    % Get number of time points
    st = size(V,1);
    st
    
    % Get repetition time
    a = V.private;
    TR = a.timing.tspace;
    TR_round = round(TR);
    
    % Only use datasets with approximate TR of 1, 2 or 3 seconds
    if (abs(TR - TR_round) < 0.05)
        
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
        
        if ( (file_number >= 924) && (file_number <= 969) )
            fMRI_volumes = permute(fMRI_volumes,[2 3 1 4]);
        end
    
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
            number_of_GPUs = 3;
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
        
        
        % Normalize
        X_GLM_1(:,1) = X_GLM_1(:,1)/norm(X_GLM_1(:,1));
        X_GLM_1(:,2) = X_GLM_1(:,2)/norm(X_GLM_1(:,2));
        
        % Orthogonalize
        X_GLM_1(:,2) = X_GLM_1(:,2) - (X_GLM_1(:,1)'*X_GLM_1(:,2))*X_GLM_1(:,1);
        
        xtxxt_GLM_1 = inv(X_GLM_1'*X_GLM_1)*X_GLM_1';
        
        c = [1 ; 0];
        
        ctxtx_GLM_1 = c'*inv(X_GLM_1'*X_GLM_1)*c;
        
        
        %----------------
        
        
        certainty = brain_voxels;
        
        smoothed_alpha_certainty = conv3(brain_voxels + eps,x_alpha_smoothing_filter);
        smoothed_alpha_certainty = conv3(smoothed_alpha_certainty,y_alpha_smoothing_filter);
        smoothed_alpha_certainty = conv3(smoothed_alpha_certainty,z_alpha_smoothing_filter);
        
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
        
        % Perform whitening on the GPU
        
        [whitened_residuals, original_residuals] = Whiten_fMRI_Volumes(fMRI_volumes, brain_voxels, certainty, f1, f2, f3, number_of_iterations, motion_compensation_filter_size, ...
            smoothing_filters_x, smoothing_filters_y, smoothing_filters_z, alpha_smoothing_filter_x, alpha_smoothing_filter_y, alpha_smoothing_filter_z, ...
            smoothing_filter_size, smoothed_certainties, smoothed_alpha_certainty, X_Detrend, xtxxt_Detrend');
        
        
        
        % Calculate spectra

        % Calculate Fourier transform of each timeseries in the mask                      
        
        TIMESERIES_ = zeros(512,1);
        for x = 1:sx
            for y = 1:sy
                for z = 1:sz
                    if (brain_voxels(y,x,z) == 1)
                        % 512 point FFT
                        timeseries = squeeze(original_residuals(y,x,z,:));
                        timeseries = timeseries / std(timeseries);
                        timeseries = fft(timeseries,512);
                        % Calculate power of spectra
                        TIMESERIES_ = TIMESERIES_ + ((abs(timeseries)).^2)/st;
                    end
                end
            end
        end
        % Calculate mean spectra
        TIMESERIES_ = TIMESERIES_ / sum(brain_voxels(:));
        
        load mean_original_spectras.mat
        mean_original_spectras(file_number,1,:) = TIMESERIES_;
        save mean_original_spectras.mat mean_original_spectras
        
        
        TIMESERIES_ = zeros(512,1);
        for x = 1:sx
            for y = 1:sy
                for z = 1:sz
                    if (brain_voxels(y,x,z) == 1)
                        % 512 point FFT
                        timeseries = squeeze(whitened_residuals(y,x,z,:));
                        timeseries = timeseries / std(timeseries);
                        timeseries = fft(timeseries,512);
                        % Calculate power of spectra
                        TIMESERIES_ = TIMESERIES_ + ((abs(timeseries)).^2)/st;
                    end
                end
            end
        end
        % Calculate mean spectra
        TIMESERIES_ = TIMESERIES_ / sum(brain_voxels(:));
        
        load mean_whitened_spectras.mat
        mean_whitened_spectras(file_number,1,:) = TIMESERIES_;
        save mean_whitened_spectras.mat mean_whitened_spectras
                        
    end
    
end


