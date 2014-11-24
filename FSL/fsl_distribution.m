close all
clear all
clc

cd /home/andek/Data/fcon1000/Beijing/

subjects = dir;

addpath('/home/andek/nifti_matlab')

all_values = [];
    
for subject = 1:198

    subjectString = subjects(subject+10).name
    
    data = load_nii(['/home/andek/Research_projects/SingleSubject/Results/Beijing/6mm/boxcar10/' subjectString '.feat/stats/tstat1.nii.gz']);
    data = data.img;
    data = double(data);
        
    [sy sx sz st] = size(data);
    
    mask = double(data ~= 0);
    brain_voxels = sum(data(:) ~= 0);
    
    all_values_temp = zeros(brain_voxels,1);
    
    v = 1;
    for z = 1:sz
        for y = 1:sy
            for x = 1:sx
                if mask(y,x,z) == 1                    
                    all_values_temp(v) = data(y,x,z);
                    v = v + 1;
                end
            end
        end
    end
    
    temp = zeros(brain_voxels + length(all_values),1);
    temp(1:brain_voxels) = all_values_temp;
    temp(brain_voxels+1:end) = all_values;      
    all_values = temp;
    size(all_values)
    
end

tvalues = trnd(217,length(all_values),1);

figure
[N1,bins] = hist(all_values,101);
N1 = N1/length(all_values);
plot(bins,N1,'r')
hold on

N2 = hist(tvalues,bins);
N2 = N2/length(tvalues);
plot(bins,N2,'g')
hold off



legend('Estimated t-values','True t-distribution')
xlabel('t-value')
ylabel('Probability')

%save('fsl_allvalues.mat','all_values');


