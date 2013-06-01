close all
clear all
clc

load ../spectra_ufp0.05_scaling_motion_regressors_no_whitening_standardized/number_of_datasets_TR
load ../spectra_ufp0.05_no_scaling_no_motion_regressors_no_whitening_standardized/rest_TRs

files = 1:1484;

load ../spectra_rpt/mean_original_spectras

% Average original residual spectra over subjects
mean_original_spectras_TR = zeros(512,3,8);

for file = files
    file
    TR = rest_TRs(file);
    round_TR = round(rest_TRs(file));
    
    % Only use datasets with integer repetition time
    if (abs(TR - round_TR) < 0.05)
        
        for experiment = 1:1
            if sum(isnan( squeeze(mean_original_spectras(file,experiment,:)) )) == 0
                mean_original_spectras_TR(:,round_TR,experiment) = mean_original_spectras_TR(:,round_TR,experiment) + squeeze(mean_original_spectras(file,experiment,:));
            end
        end
            
    else
        disp('Skipped dataset')
    end
end



load ../spectra_rpt/mean_whitened_spectras

% Average whitened residual spectra over subjects
mean_whitened_spectras_TR = zeros(512,3,8);

for file = files
    file
    TR = rest_TRs(file);
    round_TR = round(rest_TRs(file));
    
    % Only use datasets with integer repetition time
    if (abs(TR - round_TR) < 0.05)
        
        for experiment = 1:1
            if sum(isnan( squeeze(mean_whitened_spectras(file,experiment,:)) )) == 0
                mean_whitened_spectras_TR(:,round_TR,experiment) = mean_whitened_spectras_TR(:,round_TR,experiment) + squeeze(mean_whitened_spectras(file,experiment,:));
            end
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

std_whitened_spectras_TR = zeros(512,3,8);

for file = files
    file
    TR = rest_TRs(file);
    round_TR = round(rest_TRs(file));
    
    % Only use datasets with integer repetition time
    if (abs(TR - round_TR) < 0.05)
        
        for experiment = 1:1
            if sum(isnan( squeeze(mean_whitened_spectras(file,experiment,:)) )) == 0
                std_whitened_spectras_TR(:,round_TR,experiment) = std_whitened_spectras_TR(:,round_TR,experiment) + (squeeze(mean_whitened_spectras(file,experiment,:)) - squeeze(mean_whitened_spectras_TR(:,round_TR,experiment))).^2;
            end
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
    
    RPT_estimate = ORIGINAL ./ (WHITENED + eps);
    RPT_estimate(256:256+round(2*TR)) = 0;
    ratio = sum(ORIGINAL) / sum(RPT_estimate);
    RPT_estimate = ratio * RPT_estimate;
    
    f = linspace(0,(1/TR)/2,257);
    figure(TR)
    plot(f,ORIGINAL(256:end),':k')
    hold on
    plot(f,RPT_estimate(256:end),'-.k')
    hold off
    xlabel('Frequency (Hz)','FontSize',15)
    ylabel('Power','FontSize',15)
    legend('Original spectrum','Ratio between original and whitened spectra')
    if TR == 1
        title(['Average original spectrum for a repetition time of ' num2str(TR) ' second\newline        Voxel-wise AR(4) whitening prior to permutations'],'FontSize',15)
    elseif TR == 2
        title(['Average original spectrum for a repetition time of ' num2str(TR) ' seconds\newline        Voxel-wise AR(4) whitening prior to permutations'],'FontSize',15)
    elseif TR == 3
        title(['Average original spectrum for a repetition time of ' num2str(TR) ' seconds\newline        Voxel-wise AR(4) whitening prior to permutations'],'FontSize',15)
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
    
    
    filename1 = ['spectra_rpt_original_TR_' num2str(TR) '.png'];
    filename2 = ['spectra_rpt_original_TR_' num2str(TR) '.eps'];
    
    
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
    legend('Whitened spectrum','Standard deviation')
    if TR == 1 
        title(['Average whitened spectrum for a repetition time of ' num2str(TR) ' second\newline          Voxel-wise AR(4) whitening prior to permutations'],'FontSize',15)
    elseif TR == 2
        title(['Average whitened spectrum for a repetition time of ' num2str(TR) ' seconds\newline         Voxel-wise AR(4) whitening prior to permutations'],'FontSize',15)
    elseif TR == 3
        title(['Average whitened spectrum for a repetition time of ' num2str(TR) ' seconds\newline         Voxel-wise AR(4) whitening prior to permutations'],'FontSize',15)
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
    
        
    filename1 = ['spectra_rpt_whitened_TR_' num2str(TR) '.png'];
    filename2 = ['spectra_rpt_whitened_TR_' num2str(TR) '.eps'];
        
    %print('-dpng',filename1)
    %print('-depsc',filename2)

end





