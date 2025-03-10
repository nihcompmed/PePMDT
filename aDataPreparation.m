%% Load RNA-seq data from blood
clear

allData = readtable('HLLD_blood_norm_count_121s_07082018.csv');

% Column of gene names from the data
genes = table2cell(allData(:,"Symbol"));
save genes genes

%% Get individual patients data

patient1 = allData(:,contains(allData.Properties.VariableNames,"x7"));
patient2 = allData(:,contains(allData.Properties.VariableNames,"x9"));
patient3 = allData(:,contains(allData.Properties.VariableNames,"x10"));

save patient1 patient1
save patient2 patient2
save patient3 patient3

for i = 4:12
    eval(['patient',num2str(i), '= allData(:,contains(allData.Properties.VariableNames,"x',num2str(8+i),'"));'])
    eval(['save patient',num2str(i),' patient',num2str(i)])
end