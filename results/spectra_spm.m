close all
clear all
clc

load ../spectra_ufp0.05_no_scaling_no_motion_regressors_no_whitening_standardized/rest_TRs
load ../spectra_ufp0.05_no_scaling_no_motion_regressors_no_whitening_standardized/number_of_datasets_TR

files = 1:1484;

%load ../spectra_ufp0.05_scaling_motion_regressors_no_whitening_standardized/mean_spectras
load ../spectra_ufp0.05_no_scaling_no_motion_regressors_no_whitening_standardized/mean_spectras

% Average original residual spectra over subjects
mean_original_spectras_TR = zeros(512,3,8);

for file = files
    file
    TR = rest_TRs(file);
    round_TR = round(rest_TRs(file));
    
    % Only use datasets with integer repetition time
    if (abs(TR - round_TR) < 0.05)
        
        for experiment = 1:1            
            mean_original_spectras_TR(:,round_TR,experiment) = mean_original_spectras_TR(:,round_TR,experiment) + squeeze(mean_spectras(file,experiment,:));
        end
            
    else
        disp('Skipped dataset')
    end
end


%load ../spectra_ufp0.05_scaling_motion_regressors_whitening_standardized/mean_spectras
load ../spectra_ufp0.05_no_scaling_no_motion_regressors_whitening_standardized/mean_spectras
% Average whitened residual spectra over subjects
mean_whitened_spectras_TR = zeros(512,3,8);

for file = files
    file
    TR = rest_TRs(file);
    round_TR = round(rest_TRs(file));
    
    % Only use datasets with integer repetition time
    if (abs(TR - round_TR) < 0.05)
        
        for experiment = 1:1
            mean_whitened_spectras_TR(:,round_TR,experiment) = mean_whitened_spectras_TR(:,round_TR,experiment) + squeeze(mean_spectras(file,experiment,:));
        end
            
    else
        disp('Skipped dataset')
    end
end


for TR = 1:3
    for experiment = 1:1
        mean_original_spectras_TR(:,TR,experiment) = mean_original_spectras_TR(:,TR,experiment) / number_of_datasets_TR(TR,experiment);
        mean_whitened_spectras_TR(:,TR,experiment) = mean_whitened_spectras_TR(:,TR,experiment) / number_of_datasets_TR(TR,experiment);
    end
end

% Calculate standard deviation
std_whitened_spectras_TR = zeros(512,3,8);

for file = files
    file
    TR = rest_TRs(file);
    round_TR = round(rest_TRs(file));
    
    % Only use datasets with integer repetition time
    if (abs(TR - round_TR) < 0.05)
        
        for experiment = 1:1
            std_whitened_spectras_TR(:,round_TR,experiment) = std_whitened_spectras_TR(:,round_TR,experiment) + (squeeze(mean_spectras(file,experiment,:)) - squeeze(mean_whitened_spectras_TR(:,round_TR,experiment))).^2;
        end
            
    else
        disp('Skipped dataset')
    end
end


for TR = 1:3
    for experiment = 1:1
        std_whitened_spectras_TR(:,TR,experiment) = std_whitened_spectras_TR(:,TR,experiment) / number_of_datasets_TR(TR,experiment);
    end
end

for TR = 1:3
    for experiment = 1:1
        std_whitened_spectras_TR(:,TR,experiment) = sqrt(std_whitened_spectras_TR(:,TR,experiment));
    end
end


close all
experiment = 1;
for TR = 1:3
             
    ORIGINAL = squeeze(mean_original_spectras_TR(:,TR,experiment));
    ORIGINAL = fftshift(ORIGINAL);
    
    WHITENED = squeeze(mean_whitened_spectras_TR(:,TR,experiment));
    WHITENED = fftshift(WHITENED);
    
    WHITENED_STD = squeeze(std_whitened_spectras_TR(:,TR,experiment));
    WHITENED_STD = fftshift(WHITENED_STD);
    
    SPM_estimate = ORIGINAL ./ (WHITENED + eps);
    SPM_estimate(256:256+round(4*TR)) = 0;
    ratio = sum(ORIGINAL) / sum(SPM_estimate);
    SPM_estimate = ratio * SPM_estimate;
    
    f = linspace(0,(1/TR)/2,257);
    figure(TR)
    plot(f,ORIGINAL(256:end),':k')
    hold on
    plot(f,SPM_estimate(256:end),'-.k')
    hold off
    xlabel('Frequency (Hz)','FontSize',15)
    ylabel('Power','FontSize',15)
    legend('Original residuals','Ratio between original and whitened spectra')
    if TR == 1
        title(['Average original spectrum for a repetition time of ' num2str(TR) ' second\newline      SPM8, no global normalization, no motion regressors'],'FontSize',15)
        %title(['Average original spectrum for a repetition time of ' num2str(TR) ' second\newline           SPM8, global normalization, motion regressors'],'FontSize',15)        
    elseif TR == 2
        title(['Average original spectrum for a repetition time of ' num2str(TR) ' seconds\newline      SPM8, no global normalization, no motion regressors'],'FontSize',15)
        %title(['Average original spectrum for a repetition time of ' num2str(TR) ' seconds\newline           SPM8, global normalization, motion regressors'],'FontSize',15)        
    elseif TR == 3
        title(['Average original spectrum for a repetition time of ' num2str(TR) ' seconds\newline      SPM8, no global normalization, no motion regressors'],'FontSize',15)
        %title(['Average original spectrum for a repetition time of ' num2str(TR) ' seconds\newline           SPM8, global normalization, motion regressors'],'FontSize',15)        
    end
    
    if TR == 1
        set(gca,'XTick',[0:0.1:((1/TR)/2)])
    elseif TR == 2
        set(gca,'XTick',[0:0.05:((1/TR)/2)])
    else
        set(gca,'XTick',[0:0.05:((1/TR)/2)])
    end
    set(gca,'FontSize',15)
    axis([0 (1/TR)/2 0 7])    
    
    filename1 = ['spectra_no_scaling_no_motion_regressors_original_TR_' num2str(TR) '.png'];
    filename2 = ['spectra_no_scaling_no_motion_regressors_original_TR_' num2str(TR) '.eps'];
    
    %filename1 = ['spectra_scaling_motion_regressors_original_TR_' num2str(TR) '.png'];
    %filename2 = ['spectra_scaling_motion_regressors_original_TR_' num2str(TR) '.eps'];
    
    
    %print('-dpng',filename1)
    %print('-depsc',filename2)
    
    
    
    figure(TR+3)
    plot(f,WHITENED(256:end),'-k')            
    hold on
    plot(f,WHITENED(256:end)+WHITENED_STD(256:end),':k')            
    hold on
    plot(f,WHITENED(256:end)-WHITENED_STD(256:end),':k')            
    hold off
    xlabel('Frequency (Hz)','FontSize',15)
    ylabel('Power','FontSize',15)
    legend('Whitened residuals','Standard deviation')
    if TR == 1
        title(['Average whitened spectrum for a repetition time of ' num2str(TR) ' second\newline       SPM8, no global normalization, no motion regressors'],'FontSize',15)
        %title(['Average whitened spectrum for a repetition time of ' num2str(TR) ' second\newline            SPM8, global normalization, motion regressors'],'FontSize',15)        
    elseif TR == 2
        title(['Average whitened spectrum for a repetition time of ' num2str(TR) ' seconds\newline       SPM8, no global normalization, no motion regressors'],'FontSize',15)
        %title(['Average whitened spectrum for a repetition time of ' num2str(TR) ' seconds\newline            SPM8, global normalization, motion regressors'],'FontSize',15)        
    elseif TR == 3
        title(['Average whitened spectrum for a repetition time of ' num2str(TR) ' seconds\newline       SPM8, no global normalization, no motion regressors'],'FontSize',15)
        %title(['Average whitened spectrum for a repetition time of ' num2str(TR) ' seconds\newline            SPM8, global normalization, motion regressors'],'FontSize',15)        
    end
    
    if TR == 1
        set(gca,'XTick',[0:0.1:((1/TR)/2)])
    elseif TR == 2
        set(gca,'XTick',[0:0.05:((1/TR)/2)])
    else
        set(gca,'XTick',[0:0.05:((1/TR)/2)])
    end
    set(gca,'FontSize',15)
    axis([0 (1/TR)/2 0 7])    
    
    filename1 = ['spectra_no_scaling_no_motion_regressors_whitened_TR_' num2str(TR) '.png'];
    filename2 = ['spectra_no_scaling_no_motion_regressors_whitened_TR_' num2str(TR) '.eps'];
    
    %filename1 = ['spectra_scaling_motion_regressors_whitened_TR_' num2str(TR) '.png'];
    %filename2 = ['spectra_scaling_motion_regressors_whitened_TR_' num2str(TR) '.eps'];
       
    %print('-dpng',filename1)
    %print('-depsc',filename2)
end





