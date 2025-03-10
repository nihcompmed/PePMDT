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


%% Calculate Standard Deviations and Means of clusters for each time point for each patient

num_clusters = numel(final_cluster_genes_data{patient_idx});
num_timepoints = 13;
num_patients = 12;

% Initialize cell arrays to store the standard deviations and means
std_results = cell(1, num_patients); % 1 cell for each patient
mean_results = cell(1, num_patients); % 1 cell for each patient

for patient_idx = 1:num_patients
    % Initialize cells to store std and mean values for each cluster for this patient
    std_results{patient_idx} = cell(1, num_clusters);
    mean_results{patient_idx} = cell(1, num_clusters);
    
    for cluster_idx = 1:num_clusters
        % Extract gene data for this cluster (skip first column with gene names)
        cluster_data = final_cluster_genes_data{patient_idx}{cluster_idx}(:, 2:end);
        
        % Convert cell array to matrix for calculations, if necessary
        if iscell(cluster_data)
            cluster_data = cell2mat(cluster_data);
        end
        
        % Calculate standard deviation and mean across genes for each timepoint
        std_values = std(cluster_data, 0, 1);  % Std deviation along the rows (genes) for each column (timepoint)
        mean_values = mean(cluster_data, 1);   % Mean along the rows (genes) for each column (timepoint)
        
        % Store the values in the result cells for this cluster
        std_results{patient_idx}{cluster_idx} = std_values;
        mean_results{patient_idx}{cluster_idx} = mean_values;
    end
end

%% Store the standard deviation data for Python code

% Initialize an empty table to store all patients' standard deviation data
all_patients_std_data = [];

% Loop through each patient for standard deviation
for patient_idx = 1:num_patients
    % Initialize matrix for storing cluster-timepoint std for the patient
    patient_std_data = zeros(num_clusters, num_timepoints);
    
    % Loop through each cluster to fill in the matrix
    for cluster_idx = 1:num_clusters
        % Get the standard deviation values for this cluster across all timepoints
        cluster_std = std_results{patient_idx}{cluster_idx};
        patient_std_data(cluster_idx, :) = cluster_std;
    end
    
    % Convert to table and add columns for patient and cluster labels
    timepoint_labels = arrayfun(@(x) sprintf('Timepoint_%d', x), 1:num_timepoints, 'UniformOutput', false);
    T_std = array2table(patient_std_data, 'VariableNames', timepoint_labels);
    T_std.Patient = repmat({sprintf('Patient_%d', patient_idx)}, num_clusters, 1);
    T_std.Cluster = arrayfun(@(x) sprintf('Cluster_%d', x), (1:num_clusters)', 'UniformOutput', false);
    
    % Append to the main table for standard deviation
    all_patients_std_data = [all_patients_std_data; T_std]; %#ok<AGROW>
end

% Reorder columns to have Patient and Cluster as the first columns
all_patients_std_data = all_patients_std_data(:, [{'Patient', 'Cluster'}, timepoint_labels]);

% Save to CSV file
writetable(all_patients_std_data, 'input_data_std.csv');

%% Store the mean data for Python code

% Initialize an empty table to store all patients' mean data
all_patients_mean_data = [];

% Loop through each patient for mean values
for patient_idx = 1:num_patients
    % Initialize matrix for storing cluster-timepoint mean for the patient
    patient_mean_data = zeros(num_clusters, num_timepoints);
    
    % Loop through each cluster to fill in the matrix
    for cluster_idx = 1:num_clusters
        % Get the mean values for this cluster across all timepoints
        cluster_mean = mean_results{patient_idx}{cluster_idx};
        patient_mean_data(cluster_idx, :) = cluster_mean;
    end
    
    % Convert to table and add columns for patient and cluster labels
    T_mean = array2table(patient_mean_data, 'VariableNames', timepoint_labels);
    T_mean.Patient = repmat({sprintf('Patient_%d', patient_idx)}, num_clusters, 1);
    T_mean.Cluster = arrayfun(@(x) sprintf('Cluster_%d', x), (1:num_clusters)', 'UniformOutput', false);
    
    % Append to the main table for mean values
    all_patients_mean_data = [all_patients_mean_data; T_mean]; %#ok<AGROW>
end

% Reorder columns to have Patient and Cluster as the first columns
all_patients_mean_data = all_patients_mean_data(:, [{'Patient', 'Cluster'}, timepoint_labels]);

% Save to CSV file
writetable(all_patients_mean_data, 'input_data_mean.csv');

% disp('Saved all patients mean and standard deviations to input_data_mean.csv and input_data_std.csv');
