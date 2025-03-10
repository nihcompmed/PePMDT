%% Analyse initial data

% 1. First observe the existing data timepoint wise for each patient

% 2. Select patients data for training FNN based on the amount of data
% present

% 3. Fill the missing values

%% Data Time Points

% Patients' time points(in minutes): PreOp | IntraOp_Pre_resection |
% | 5 min | 30 min | 60 min | 120 min |
% | 180 min (3hrs)| 240 min (4hrs)| 1440 min (1day)| 2880 min (2days)|
% | 4320 min (3days)| 5760 min (4days)| 14400 min (10days)| 129600 min
% (3months)| 259200 min (6months)| 525600 min (1year)|

% full names for Day1 and Day10 are written to avoid mistake
TP = {'01PreOp','Pre_resection','5min','30min','60min','120min','07PostOp_Day1',...
    'Day2','Day3','Day4','12PostOp_Day10','Month3','Month6','Year1'};
% No patient is having 180min and 240min time points but they are mentioned
% in the article, that is why were not considered for the data

%% Load initial DEG data (2-fold change at any timepoint wrt initial preOp data)

% load patients data
num_patients = 12;

for i = 1:num_patients
    load(['DEGpatient',num2str(i),'.mat'])
end

%% Find missing time points 

xq = 1:1:14; % time points

patientstimepoints = (1:num_patients)';
for i = 1:num_patients
    eval(['rawData = DEGpatient',num2str(i),';']) % raw data
    data = table2array(rawData); % original data

    for j = 1:length(TP)
        eval(['tf = any(contains(DEGpatient',num2str(i),'.Properties.VariableNames,TP{j}));'])
        if double(tf) == 1
            patientstimepoints(i,j+1) = xq(j);
        end
    end
end

%% Sort the patients according to the data available to us (most timepoints --> least)

% Extract the relevant portion of the matrix (excluding the first column)
data = patientstimepoints(:, 2:end);

% Step 1: Count the non-zero elements in each row
non_zero_counts = sum(data ~= 0, 2);

% Step 2: Sort the rows based on the non-zero element count in decreasing order
[~, sorted_indices] = sort(non_zero_counts, 'descend');

% Step 3: Re-arrange the matrix based on the sorted indices
sorted_data = data(sorted_indices, :);

% Optionally, re-arrange the entire patientstimepoints matrix, including the first column
sorted_patientstimepoints = patientstimepoints(sorted_indices, :);

%% Plot

% Create a binary matrix with 0s and 1s based on non-zero elements
binary_matrix = sorted_patientstimepoints(:, 2:end) ~= 0;

% Plot the binary matrix
figure;
imagesc(binary_matrix);
colormap([1 1 1; 0 0.5 0.5]); % Set the colormap: white for 0, teal for 1
xlabel('Timepoints');
ylabel('Patients');

ax = gca;
ax.FontSize = 26;
ax.LineWidth = 1.8;

% Set X and Y ticks for labels (use integers for centering)
ax.XTick = 1:size(binary_matrix, 2);  % Integer ticks for timepoints
ax.YTick = 1:12;  % Integer ticks for patients

% Set custom x-tick labels for timepoints (centered on the integer ticks)
ax.XTickLabel = {'Pre Operation','Before surgery','5 minutes','30 minutes','60 minutes','120 minutes',...
    '1 day','2 days','3 days','4 days','10 days','3 months','6 months','1 year'};

% Set custom y-tick labels for patients (centered on the integer ticks)
for i = 1:12
    ax.YTickLabel{i} = ['Patient ' num2str(sorted_patientstimepoints(i,1))];
end

% Adjust grid to outline the boxes, placed between the tick marks
set(gca, 'XTick', 0.5:1:(size(binary_matrix, 2)+0.5), 'YTick', 0.5:1:12.5, 'TickLength', [0 0]);

% Enable the grid and set grid properties (grid around the boxes)
grid on;
ax.GridColor = [0 0 0];  % Black grid lines
ax.GridAlpha = 1;        % Fully visible grid lines
ax.GridLineStyle = '-';  % Solid lines

% Ensure that tick marks do not interfere with the visual appearance
set(gca, 'TickLength', [0 0]);

% Ensure the tick labels remain on top of the grid
set(gca, 'Layer', 'top');

% --- Add counts of 1s at the right-hand side ---
for i = 1:12
    % Count the number of 1s for each patient (row)
    count_of_ones = sum(binary_matrix(i, :));
    total_timepoints = size(binary_matrix, 2);  % Total timepoints

    % Display the count at the right side of the plot (slightly outside the matrix)
    text(total_timepoints + 0.5, i, ['  (',num2str(count_of_ones) '/' num2str(total_timepoints),')'], ...
        'FontSize', 20, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
end


%% Now fill the missing values with the average values
% so this will fill the missing values based on the binary
% matrix/sorted_patienttimepoints matrix


% First fill the missing columns with zeros
for i = 1:num_patients
    eval(['rawData = DEGpatient',num2str(i),';']) % raw data
    data = table2array(rawData); % original data

    % Create a new matrix for each patient
    eval(['filledDEGpatient',num2str(i),' = zeros(size(data,1),14);']) % Initialize with zeros

    col_no = 0;
    for j = 2:size(patientstimepoints,2)
        if find(patientstimepoints(i,j)) == 1
            col_no = col_no+1;
            eval(['filledDEGpatient',num2str(i),'(:,patientstimepoints(i,j))=data(:,col_no);'])
        end
    end    
end

% Now based on the binary/sorted_patienttimepoints matrix fill the missing columns with average values

% Get the sum of the data points for each of the training data matrices
% Initialize the accumulator matrix with zeros (assuming the size of the matrices is the same)
sumData = zeros(size(filledDEGpatient1));


%% This is for all datasets : fill missing values for all patients

for i = 1:num_patients
    % Get the patient number from the sorted list
    patient_num = sorted_patientstimepoints(i, 1);
    
    % Access the corresponding filledDEGpatient matrix
    eval(['patient_data = filledDEGpatient', num2str(patient_num), ';']);
    
    % Sum the matrices
    sumData = sumData + patient_data;
end


for i = 1:num_patients
    for j = 2:size(sorted_patientstimepoints,2)
        if sorted_patientstimepoints(i,j) == 0
            [indx,~] = find(sorted_patientstimepoints(1:8,j));
            meanData = sumData(:,j-1) / length(indx);
            eval(['filledDEGpatient',num2str(sorted_patientstimepoints(i,1)),'(:,j-1)=meanData;'])
        end
    end
end

for i = 1:12
    eval(['save(''filledDEGpatient', num2str(i), '.mat'', ''filledDEGpatient', num2str(i), ''')']); 
end

