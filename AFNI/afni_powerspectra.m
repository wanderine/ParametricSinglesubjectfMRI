close all
clear all
clc

cd /home/andek/Data/fcon1000/Beijing/

subjects = dir;

addpath('/home/andek/nifti_matlab')

MEAN_TIMESERIES = zeros(512,1);
    
for subject = 1:198

    subjectString = subjects(subject+10).name
    
    data = load_nii(['/home/andek/Research_projects/SingleSubject/Results/Beijing/6mm/boxcar10_REML/' subjectString '.results/whitened_residuals.nii.gz']);
    data = data.img;
    data = double(data);
    
    mask = load_nii(['/home/andek/Research_projects/SingleSubject/Results/Beijing/6mm/boxcar10_REML/' subjectString '.results/mask.nii.gz']);
    mask = mask.img;
    mask = double(mask);

    [sy sx sz st] = size(data);
    
    brain_voxels = sum(mask(:));
    MEAN_TIMESERIES_TEMP = zeros(512,1);
    
    for z = 1:sz
        for y = 1:sy
            for x = 1:sx
                if mask(y,x,z) == 1                    
                    timeseries = squeeze(data(y,x,z,:));
                    if (std(timeseries) ~= 0)
                        timeseries = timeseries/(std(timeseries) + eps);
                        TIMESERIES = fft(timeseries,512);
                        MEAN_TIMESERIES_TEMP = MEAN_TIMESERIES_TEMP + ((abs(TIMESERIES)).^2)/st;
                    end
                end
            end
        end
    end
      
    MEAN_TIMESERIES_TEMP = MEAN_TIMESERIES_TEMP / brain_voxels;
    MEAN_TIMESERIES = MEAN_TIMESERIES + MEAN_TIMESERIES_TEMP;
    
end

MEAN_TIMESERIES = MEAN_TIMESERIES / 198;
MEAN_TIMESERIES = fftshift(MEAN_TIMESERIES);



%save('afni_powerspectra.mat','MEAN_TIMESERIES');

