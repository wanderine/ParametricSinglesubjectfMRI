
clear all
clc
close all

addpath('D:\spm8')
addpath('C:\Users\wande\Documents\GitHub\ParametricSinglesubjectfMRI\SPM')

study = 'Cambridge'

cd D:\fcon1000\Cambridge\
subjects = dir;

smoothnessEstimates = zeros(198,1);
    
analyses = 0;

for subject = 1:198

    analyses = analyses + 1;   
    subjectString = subjects(subject+2).name;
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASSICAL STATISTICAL ANALYSIS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    for experiment = 1
        for smoothing = 6
            
            subjectString
                                              
            if experiment == 1
                % Boxcar 10
                experimentString = 'boxcar10';
            elseif experiment == 2
                % Boxcar 30
                experimentString = 'boxcar30';
            end
                                                            
            % Load SPM file
            clear SPM
            load(['D:\fcon1000\' study '\' subjectString '\func\SPM_' experimentString '_s' num2str(smoothing) '.mat']);
            
            FWHM = SPM.xVol.FWHM;
            M    = SPM.xVol.M;
            FWHM(1) = FWHM(1) * abs(M(1,1));
            FWHM(2) = FWHM(2) * abs(M(2,2));
            FWHM(3) = FWHM(3) * abs(M(3,3));
            smoothnessEstimates(subject) = FWHM(1);
                        
        end
    end   
end

