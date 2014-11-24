close all
clear all
clc

cd D:\fcon1000\Beijing

study = 'Beijing';
experimentString = 'boxcar10';
smoothing = 6;

subjects = dir;

all_values = [];
    
for subject = 1:198

    subjectString = subjects(subject+2).name
    
    V = spm_vol(['D:\fcon1000\' study '\' subjectString '\func\spmT_' experimentString '_s' num2str(smoothing) '.hdr']);
    [data,aa] = spm_read_vols(V);
                
    mask = double(data ~= 0);
    brain_voxels = sum(data(:) ~= 0);
    
    all_values_temp = zeros(brain_voxels,1);
    
    [sy sx sz st] = size(data);
    
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

%save('spm_allvalues.mat','all_values');


