%% WGCNA with average data

load DEGgenes.mat

sumData = zeros(length(DEGgenes),14);

% Loop over each patient (load data for each patient and process it)
for patient = setdiff(1:12,[5,6,10,11])
    % Load the corresponding patient data file (e.g., 'DEGpatient1.mat', 'DEGpatient2.mat', ...)
    filename = ['filledDEGpatient', num2str(patient), '.mat'];
    load(filename);  % This loads the variable DEGpatientX into the workspace

    % Extract the data for fold change calculation (columns 3:end for each patient)
    patient_data = eval(['filledDEGpatient', num2str(patient)]);
    % Sum the matrices
    sumData = sumData + patient_data;
end

meanData = sumData/length(setdiff(1:12,[5,6,10,11]));

data = meanData;

expression_data = data(:, 3:end);  % Extract all columns except the 1st (which is the baseline)
baseline = data(:, 1);             % Column 1 is the baseline for fold changes

% Step 1: Calculate fold changes for columns 3 to end with respect to the baseline
fold_changes = expression_data ./ baseline;  % Element-wise division for fold change calculation

% Step 2: Identify 2-fold up or down regulation for two consecutive columns
two_fold_up = fold_changes >= 2;     % Logical array where fold change >= 2 (upregulation)
two_fold_down = fold_changes <= 0.5; % Logical array where fold change <= 0.5 (downregulation)

% Step 3: Find rows with 2-fold change (up or down) for two consecutive time points
consecutive_up = two_fold_up(:, 1:end-1) & two_fold_up(:, 2:end);  % Check consecutive columns
consecutive_down = two_fold_down(:, 1:end-1) & two_fold_down(:, 2:end);

% Step 4: Combine the up and down conditions
consecutive_changes = consecutive_up | consecutive_down;

% Step 5: Find the row indices of genes that satisfy this condition for the current patient
[row_indices, ~] = find(consecutive_changes);

% Step: Remove rows where the fold change between the baseline and the second column
% is greater than 1.5-fold up or down
fold_change_second_column = data(:, 2)./ baseline;  % Fold changes for the second column
rows_to_remove = fold_change_second_column > 1.5 | fold_change_second_column < 0.67;  % Rows with >1.5-fold change

% Store the unique row indices for this patient
degs_allpatients = setdiff(unique(row_indices),find(rows_to_remove));

%% compare DEG genes patientwise with Lawrence et al paper results

% Load the data from the Excel file
filename = 'LawrencePaperGeneData.xlsx';
sheets = {'Peak 1', 'Peak 2', 'Peak 3'}; % Names of the sheets
excluded_patients = [5,6,10,11]; % Excluded patients
patients = setdiff(1:12, excluded_patients); % Remaining patients for comparison

% Initialize a structure to hold the results for each patient
patientResults = struct();

% Get the DEG genes for the current patient
deggenes = DEGgenes(degs_allpatients);
DEGgenes = DEGgenes(degs_allpatients);

% Initialize a results cell array for this patient
results = {};

% Loop through each peak (sheet)
for sheetIdx = 1:length(sheets)
    sheetName = sheets{sheetIdx};

    % Read the data from the current sheet
    [~, ~, rawData] = xlsread(filename, sheetName);

    % Loop through each row of the sheet (starting from row 1)
    for rowIdx = 1:size(rawData, 1)
        % The first column contains the function name
        functionName = rawData{rowIdx, 1}; % Function name

        % Get the gene names in the current row, excluding the first column
        associatedGenes = rawData(rowIdx, 2:end); % Genes for the current function

        % Remove any NaN or empty cells
        associatedGenes = associatedGenes(~cellfun(@(x) isempty(x) || (isnumeric(x) && isnan(x)), associatedGenes));

        % Find common genes between the patient's deggenes and the associated genes
        commonGenes = intersect(deggenes, associatedGenes);

        if ~isempty(commonGenes) % If there are common genes
            % Store the results for the current patient
            results = [results; {sheetName, functionName, strjoin(commonGenes, ', ')}];
        end
    end
end

fprintf('Compare all training patients DEGs with Lawrence Paper: \n\n')
resultTable = cell2table(results, ...
    'VariableNames', {'Peak', 'Function', 'Common_Genes'});
disp(resultTable);

%% WGCNA for average patient data

noRows = 4;
noColumns = 4;

% Extract the data (excluding baseline column)
data = meanData(:, 2:end)./meanData(:, 1);

% Take the log2 of the resulting data
log2_data = log2(data + 1);

%% Normalization

% data = log2_data;
%
% % Initialize a matrix to store the normalized data
% zscore_normalized_data = zeros(size(data));
%
% % Apply z-score normalization row-wise
% for i = 1:size(data, 1)
%     row_data = data(i, :);  % Extract the i-th row
%     row_mean = mean(row_data);  % Compute the mean of the row
%     row_std = std(row_data);    % Compute the standard deviation of the row
%
%     % Perform z-score normalization for the i-th row
%     zscore_normalized_data(i, :) = (row_data - row_mean) / row_std;
% end
%
% data = zscore_normalized_data(deg_pat12,:);

%% Perform PCA

data = log2_data(degs_allpatients,:);

[pc, scores, pcvars] = pca(data);

% Compute the exact percentage of the variance accounted for by each component
pc_per = pcvars./sum(pcvars) * 100;

% cumsum command to see the cumulative sum of the variances
pc_cum = cumsum(pcvars./sum(pcvars) * 100);

%% PCA Figure

figure
scatter(scores(:,1),scores(:,2), 'LineWidth', 2.5);

ax = gca;
ax.FontSize = 24;
ax.LineWidth = 1.8;

% ax.XLim = [-250 150];
% ax.XTick= -250:50:150;
% ax.XLim = [-150 150];
% ax.XTick= -150:50:150;
% ax.YLim = [-150 150];
% ax.YTick= -150:50:150;

xlabel('PC1');
ylabel('PC2');
title('Principal Component Scatter Plot');

%% Self-Organizing Maps (SOM) using Deep Learning Toolbox

% The selforgmap function creates a new SOM network object.
% generate a SOM using the first two principal components.

P = scores(:,1:2)';
net = selforgmap([noRows noColumns]);

% rng(42);                    % Fix the seed

% Train the network using the default parameters.

net = train(net,P);

% Use plotsom to display the network over a scatter plot of the data.
% SOM algorithm uses random starting points so the results will vary from run to run.

%% Figure

figure;
plot(P(1,:),P(2,:),'.g','markersize',30)
hold on
plotsom(net.iw{1,1},net.layers{1}.distances)
hold off

ax = gca;
ax.FontSize = 24;
ax.LineWidth = 1.8;

% We can assign clusters using the SOM by finding the nearest node to each
% point in the data set.

distances = dist(P',net.IW{1}');
[d,cndx] = min(distances,[],2); % cndx contains the cluster index

fig1 = figure('Visible', 'on', 'Position', [0, 0, 2000, 1000]); % Adjust the size as needed

gscatter(P(1,:),P(2,:),cndx,hsv(numel(unique(cndx))),"o", "filled");
% legend off;
hold on
plotsom(net.iw{1,1},net.layers{1}.distances);
hold off

ax = gca;
ax.FontSize = 24;
ax.LineWidth = 1.8;

%% Gene Trajectories of each cluster: Mean and std

cluster_data_mean = [];

clr = hsv(numel(unique(cndx)));

fig2 = figure('Visible', 'on', 'Position', [0, 0, 2000, 1000]); % Adjust the size as needed

for i = 1:numel(unique(cndx))

    % Mean and std of multi-patient gene clusters
    mm = mean(data(cndx==i,:),1);
    cluster_data_mean = [cluster_data_mean;mm];
    sd = std(data(cndx==i,:));

    % Plot clusters
    subplot(noRows,noColumns,i)
    fill([1:length(mm) fliplr(1:length(mm))], [mm+sd fliplr(mm-sd)],...
        [eval(strcat('clr(',num2str(i),',:)'))],'linestyle', 'none','FaceAlpha',0.5)
    hold on
    plot(1:length(mm),mm,'--ok','MarkerFaceColor','k','markersize',8,'LineWidth',2.5)

    ax = gca;
    ax.FontSize = 16;
    ax.LineWidth = 1.8;
    ax.XLim = [0.5 length(mm)+0.5];
    ax.XTick= 1:1:length(mm);

    ax.XTickLabel = {'Before resection','5 minutes','30 minutes', '60 minutes','120 minutes',...
        '1 day','2 days','3 days','4 days','10 days','3 months','6 months','1 year'};

    title(strcat('Cluster',{' '},num2str(i)));

    % set(gca,'FontWeight','bold', 'FontName','Calibri')
    set(gca,'TickLength',[0.005, 0.01])
end

%% Cluster data

for i = 1:numel(unique(cndx))
    Cluster_dataAllPat{i} = [DEGgenes(cndx==i), num2cell(meanData(cndx==i,:))];
end

save Cluster_dataAllPat Cluster_dataAllPat

%% Match Lawrence paaper results with cluster data

% Load the data from the Excel file
filename = 'LawrencePaperGeneData.xlsx'; % Ensure the correct filename and extension
sheets = {'Peak 1', 'Peak 2', 'Peak 3'}; % Names of the sheets
numClusters = noRows*noColumns; % Total number of clusters

% Initialize a cell array to hold the results
results = {}; % Use a cell array for easy appending

for clusterIdx = 1:numClusters
    % Get the list of genes in the current cluster
    clusterGenes = Cluster_dataAllPat{clusterIdx}(:, 1); % First column contains gene names

    % Loop through each peak (sheet)
    for sheetIdx = 1:length(sheets)
        sheetName = sheets{sheetIdx};

        % Read the data from the current sheet
        [~, ~, rawData] = xlsread(filename, sheetName);

        % Loop through each row of the sheet (starting from row 1)
        for rowIdx = 1:size(rawData, 1)
            % The first column contains the function name
            functionName = rawData{rowIdx, 1}; % Function name

            % Get the gene names in the current row, excluding the first column
            associatedGenes = rawData(rowIdx, 2:end); % Genes for the current function

            % Remove any NaN or empty cells
            associatedGenes = associatedGenes(~cellfun(@(x) isempty(x) || (isnumeric(x) && isnan(x)), associatedGenes));

            % Find common genes between the cluster and associated genes
            commonGenes = intersect(clusterGenes, associatedGenes);

            if ~isempty(commonGenes) % If there are common genes
                % Store the results
                results = [results; {clusterIdx, sheetName, functionName, strjoin(commonGenes, ', ')}];
            end
        end
    end
end

% Convert the results cell array to a table
resultTable = cell2table(results, 'VariableNames', {'Cluster', 'Peak', 'Function', 'Genes'});

% Display the results table
fprintf('Compare all cluster data of training patients DEGs with Lawrence Paper: \n\n')
disp(resultTable);

%% Save the WGCNA gene clusters in an excell file

% Define the name of the Excel file to write to
excelFileName = 'cluster_gene_tarining_patients.xlsx';

% Loop over each cluster in the cell array
for i = 1:length(Cluster_dataAllPat)
    % Get the gene names for the current cluster
    currentClusterGenes = Cluster_dataAllPat{i}(:,1);

    % Create a sheet name (e.g., 'Cluster 1', 'Cluster 2', etc.)
    sheet_name = ['Cluster ' num2str(i)];

    % Write the gene names to the Excel file
    writetable(cell2table(currentClusterGenes), excelFileName, 'Sheet', sheet_name, 'WriteVariableNames', false);
end


%% Cluster data for training patients | Data used ---> log2(fold-change wrt baseline) same used as clusters

% Define the training patients
training_patients = setdiff(1:12,[5,6,10,11]);

% Loop through the training patients
for i = 1:numel(training_patients)
    patient_idx = training_patients(i);  % Get the actual patient index
    eval(['patient_data = filledDEGpatient', num2str(patient_idx), '(degs_allpatients,:);'])
    fold_change_data = patient_data(:,2:end)./patient_data(:,1);
    log2_data = log2(fold_change_data+1);

    % Loop through the unique clusters
    for j = 1:numel(unique(cndx))
        % Collect the DEG genes and their corresponding data for this cluster
        eval(['cluster_dataTrainingPatients',num2str(patient_idx),'{j} = [DEGgenes(cndx==j), num2cell(log2_data(cndx==j,:))];'])
    end

    % Save the data for the current patient
    eval(['save cluster_dataTrainingPatients',num2str(patient_idx),' cluster_dataTrainingPatients',num2str(patient_idx)])
end

%% Cluster data for test patients | Data used ---> log2(fold-change wrt baseline) same used as clusters
% the missing entries in test patients will be 0 after log2(data+1) also

test_patients = [5, 6, 10, 11];

% Loop through the test patients
for i = 1:numel(test_patients)
    patient_idx = test_patients(i);  % Get the actual patient index

    filename = ['filledDEGpatient', num2str(patient_idx), '.mat'];
    load(filename);

    eval(['patient_data = filledDEGpatient', num2str(patient_idx), '(degs_allpatients,:);'])
    
    % Calculate log2 fold change
    fold_change_data = patient_data(:, 2:end) ./ patient_data(:, 1);
    log2_data = log2(fold_change_data + 1);  % +1 for stability in case of 0 values

    % Loop through the unique clusters
    for j = 1:numel(unique(cndx))
        % Collect the DEG genes and their corresponding data for this cluster
        eval(['cluster_dataTestPatients', num2str(patient_idx), '{j} = [DEGgenes(cndx==j), num2cell(log2_data(cndx==j,:))];'])
    end

    % Save the data for the current test patient
    eval(['save cluster_dataTestPatients', num2str(patient_idx), ' cluster_dataTestPatients', num2str(patient_idx)])
end


% Save the whole workspace
% save WGCNA_workspace


