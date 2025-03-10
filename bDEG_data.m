%% Differentially expressed geenes across patients
% A gene is considered as DEG if the gene is expressed 2-fold up or down
% wrt the pre operation baseline value at least for one time point during 
% post surgical recovery

% This DEGs are simply 2-folds up or down genes that are for any time
% points wrt the initial baseline timepoint

clear

%% Load gene data

load genes.mat % load gene list

% Find indices of non-empty cells
nonEmptyIdx = find(~cellfun('isempty', genes));
nonEmptyGenes = genes(nonEmptyIdx);

% Find unique cells among non-empty cells and their original indices
[uniqueGenes, uniqueIdxWithinNonEmpty] = unique(nonEmptyGenes, 'stable');
uniqueIdx = nonEmptyIdx(uniqueIdxWithinNonEmpty);

%% Load patients data

num_patients = 12;

for i = 1:num_patients
    load(['patient',num2str(i),'.mat'])
    eval(['patient',num2str(i),' = patient',num2str(i),'(uniqueIdx,:);'])
end

%% Differentially expressed genes

data = [];
for i = 1:num_patients
    eval(['tmp = table2array(patient',num2str(i),');'])
    data = [data,tmp(:,2:end)./tmp(:,1)]; % calculate fold change
end

% Filtering the genes/data

% removing NaN and inf values
nanIndices = find(any(isnan(data),2)); % get NaN row indices
infIndices = find(any(isinf(data),2)); % get Inf row indices

% find the rows with no fold change greater than 2 from that subtract the
% rows with fold change less than 0.5 by set difference
nonDEGIndices = setdiff(find(~any(data>=2,2)),find(any(data<=1/2,2))); % get nonDEG indices

naninfIndices = unique([nanIndices;infIndices;nonDEGIndices]);

% check for DEG
data(naninfIndices,:) = [];

%% Filter pateint data and gene: Consider only DEGs

uniqueGenes(naninfIndices) = [];
% For these indices data is missing for a patricular time point for every
% patient, so we need to remove these genes also. This is found when
% normalising the data for each time point.
missingRows = [3336,7111,5211,8531];
uniqueGenes(missingRows) = [];
DEGgenes = uniqueGenes;
save DEGgenes DEGgenes

for i = 1:num_patients
    eval(['patient',num2str(i),'(naninfIndices,:) = [];'])
    eval(['patient',num2str(i),'(missingRows,:) = [];'])
    eval(['DEGpatient',num2str(i),' = patient',num2str(i),';'])
    eval(['save DEGpatient',num2str(i),' DEGpatient',num2str(i)])
end

