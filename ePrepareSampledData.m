%% This code is for preparing sampled data for VAE training

% Data preparation for python codes
% Here we will consider data of 8 patients that are chosen for training

% Load data
load simulated_data.mat

% Patient numbers
training_patients = setdiff(1:12,[5,6,10,11]);
test_patients = [5,6,10,11];

for i = 1:numel(training_patients)
    patient_idx = training_patients(i);  % Get the actual patient index
    eval(['load cluster_dataTrainingPatients',num2str(patient_idx)])
    eval(['final_cluster_genes_data{patient_idx} = cluster_dataTrainingPatients',num2str(patient_idx),'([1:10, 12:end]);'])
end

for i = 1:numel(test_patients)
    patient_idx = test_patients(i);  % Get the actual patient index
    eval(['load cluster_dataTestPatients',num2str(patient_idx)])
    eval(['final_cluster_genes_data{patient_idx} = cluster_dataTestPatients',num2str(patient_idx),'([1:10, 12:end]);'])
end


%% Data for training dataset

num_clusters = length(cluster_dataTrainingPatients1)-1; % omitted cluster 11
num_samples = 8000;

samples = samplePatients(training_patients, num_samples);

[input_data_random, sim_data_random] = prepareLargerMatrix(final_cluster_genes_data, samples, num_clusters, num_samples, simulated_data);

% Generate filenames
input_filename = sprintf('train_input_data.csv');
output_filename = sprintf('train_output_data.csv');
patient_ids_filename = sprintf('train_patient_ids.csv');

% Save data to CSV files
writematrix(input_data_random, input_filename);
writematrix(sim_data_random, output_filename);
writematrix(samples, patient_ids_filename);

%% Data for test dataset

num_clusters2 = length(cluster_dataTrainingPatients1)-1;
num_samples2 = 4000;

samples2 = samplePatients(test_patients, num_samples2);

[input_data_rest, sim_data_rest] = prepareLargerMatrix(final_cluster_genes_data, samples2, num_clusters2, num_samples2, simulated_data);

% Generate filenames
input_filename = sprintf('test_input_data.csv');
output_filename = sprintf('test_output_data.csv');
patient_ids_filename = sprintf('test_patient_ids.csv');

% Save data to CSV files
writematrix(input_data_rest, input_filename);
writematrix(sim_data_rest, output_filename);
writematrix(samples2, patient_ids_filename);

