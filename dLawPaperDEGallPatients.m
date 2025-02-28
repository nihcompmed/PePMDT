% Initialize an empty cell array to store indices of genes from all patients
all_patient_degs = cell(1, 12);
load DEGgenes.mat

% Loop over each patient (load data for each patient and process it)
for patient = 1:12 %setdiff(1:12,[5,6,10,11])
    % Load the corresponding patient data file (e.g., 'DEGpatient1.mat', 'DEGpatient2.mat', ...)
    filename = ['filledDEGpatient', num2str(patient), '.mat'];
    load(filename);  % This loads the variable DEGpatientX into the workspace

    % Extract the data for fold change calculation (columns 3:end for each patient)
    data = eval(['filledDEGpatient', num2str(patient)]);  % Adjust variable based on your data structure
    expression_data = data(:, 3:end);  % Extract all columns except the 1st (which is the baseline)
    baseline = data(:, 1);             % Column 1 is the baseline for fold changes

    % Step 1: Calculate fold changes for columns 3 to end with respect to the baseline
    fold_changes = expression_data ./ baseline;  % Element-wise division for fold change calculation

    % Step 2: Identify 2-fold up or down regulation for two consecutive columns
    two_fold_up = fold_changes >= 2;     % Logical array where fold change >= 2 (upregulation)
    two_fold_down = fold_changes <= 0.5; % Logical array where fold change <= 0.5 (downregulation)

    % Step 3: Find rows with 2-fold change (up or down) for two consecutive time points
    consecutive_up = two_fold_up(:, 2:end-1) & two_fold_up(:, 3:end);  % Check consecutive columns
    consecutive_down = two_fold_down(:, 2:end-1) & two_fold_down(:, 3:end);

    % Step 4: Combine the up and down conditions
    consecutive_changes = consecutive_up | consecutive_down;

    % Step 5: Find the row indices of genes that satisfy this condition for the current patient
    [row_indices, ~] = find(consecutive_changes);

    % Step: Remove rows where the fold change between the baseline and the second column
    % is greater than 1.5-fold up or down
    fold_change_second_column = data(:, 2)./ baseline;  % Fold changes for the second column
    rows_to_remove = fold_change_second_column > 1.5 | fold_change_second_column < 0.67;  % Rows with >1.5-fold change

    % Store the unique row indices for this patient
    all_patient_degs{patient} = setdiff(unique(row_indices),find(rows_to_remove));

end

% Step 6: Find common genes across training data patients
common_genes = all_patient_degs{1};  % Start with the first patient's genes
for patient = setdiff(2:12,[5,6,10,11])
    common_genes = intersect(common_genes, all_patient_degs{patient});  % Find common genes
end

save all_patient_degs all_patient_degs

%% compare DEG genes patientwise with Lawrence et al paper results

for patient = setdiff(1:12,[5,6,10,11])
    deggenes = DEGgenes(all_patient_degs{patient});
end

% Load the data from the Excel file
filename = 'LawrencePaperGeneData.xlsx'; 
sheets = {'Peak 1', 'Peak 2', 'Peak 3'}; % Names of the sheets
excluded_patients = [5,6,10,11]; % Excluded patients
patients = setdiff(1:12, excluded_patients); % Remaining patients for comparison

% Initialize a structure to hold the results for each patient
patientResults = struct();

for patient = patients
    % Get the DEG genes for the current patient
    deggenes = DEGgenes(all_patient_degs{patient});
    
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
                results = [results; {patient, sheetName, functionName, strjoin(commonGenes, ', ')}];
            end
        end
    end
    
    % Save the results for the current patient in the structure
    patientResults(patient).results = results;
end

% Display results for each patient
for patient = patients
    fprintf('\nResults for Patient %d:\n', patient);
    
    % Convert the results to a table for easier viewing
    if ~isempty(patientResults(patient).results)
        resultTable = cell2table(patientResults(patient).results, ...
            'VariableNames', {'Patient', 'Peak', 'Function', 'Common_Genes'});
        disp(resultTable);
    else
        disp('No common genes found.');
    end
end


