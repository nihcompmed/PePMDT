% Define baseline parameter values
params_baseline = struct();
params_baseline.M = 5.8; % human donor, optimized
params_baseline.ktnf = 1.5; % TNF production
params_baseline.kaptnf = 0.9; % TNF decay
params_baseline.vjak = 2e4; % JAK production
params_baseline.kmjak = 1e4;
params_baseline.kapjak = 4e-1; % JAK decay
params_baseline.kprost3 = 2; % proSTAT3
params_baseline.vst3 = 750; % STAT3 production
params_baseline.kmst3 = 0.4;
params_baseline.kapst3 = 0.1; % STAT3 decay
params_baseline.vsoc = 2.4e4; % SOCS3 production
params_baseline.kmsoc = 7e-4;
params_baseline.kisoc = 15e-3; % SOCS3 inhibition constant
params_baseline.kapsoc = 0.4; % SOCS3 decay
params_baseline.vie = 2.5e2; % IE production
params_baseline.kmie = 18;
params_baseline.kapie = 5; % IE decay
params_baseline.kup = 6e-2; % GF uptake by ECM
params_baseline.kgf = 1.125e-1; % GF production
params_baseline.kapgf = 0.23; % GF decay
params_baseline.kapecm = 33; % ECM decay
params_baseline.kdeg = 7; % ECM degradation by TNF
params_baseline.kecm = params_baseline.kapecm + params_baseline.kdeg; % ECM production
params_baseline.kq = 7e-3; % Q -> P
params_baseline.kr = 5.4e-2; % R -> Q
params_baseline.kp = 4.4e-3; % P -> R
params_baseline.kprol = 2e-2; % Proliferation rate
params_baseline.kreq = 1e-1; % Requiescence rate
params_baseline.thetar = 8;
params_baseline.betar = 3;
params_baseline.kapop = 1e-1 / 24; % Apoptosis
params_baseline.thetaa = 0.15;
params_baseline.betaa = 0.075;

% Parameter bounds (±20%)
param_names = fieldnames(params_baseline);
lb = zeros(1, numel(param_names));
ub = zeros(1, numel(param_names));
for i = 1:numel(param_names)
    val = params_baseline.(param_names{i});
    lb(i) = val * 0.8;
    ub(i) = val * 1.2;
end

% Sampling and simulation settings
numSamples = 1000;
lhs_matrix = lhsdesign(numSamples, numel(param_names));
param_samples = repmat(lb, numSamples, 1) + lhs_matrix .* (repmat(ub - lb, numSamples, 1));
y0 = [1, 1, 1, 1, 1, 1, 1, 0.35, 0, 0];
time = linspace(0, 365*24, 525601);

% Timepoint values and names
time_points = [1,5,30,60,120,1440,2880,4320,5760,14400,129600,259200,525600];
timepoint_names = { 'Before Surgery', '5 Minutes', '30 Minutes', '60 Minutes', ...
    '120 Minutes', '1 Day', '2 Days', '3 Days', ...
    '4 Days', '10 Days', '3 Months', '6 Months', '1 Year' };

% Labels
param_names_latex = {'\it{M}', '\it{k_{TNF}}', '\it{\kappa_{TNF}}', '\it{V_{JAK}}', '\it{k_{M}^{JAK}}', ...
 '\it{\kappa_{JAK}}', '\it{[proSTAT3]}', '\it{V_{ST3}}', '\it{k_M^{ST3}}', '\it{\kappa_{ST3}}', ...
 '\it{V_{SOCS3}}', '\it{k_M^{SOCS3}}', '\it{k_I^{SOCS3}}', '\it{\kappa_{SOCS3}}', '\it{V_{IE}}', ...
 '\it{k_M^{IE}}', '\it{\kappa_{IE}}', '\it{k_{up}}', '\it{k_{GF}}', '\it{\kappa_{GF}}', ...
 '\it{\kappa_{ECM}}', '\it{k_{deg}}', '\it{k_{ECM}}', '\it{k_Q}', '\it{k_R}', ...
 '\it{k_P}', '\it{k_{prol}}', '\it{k_{req}}', '\it{\theta_{req}}', '\it{\beta_{req}}', ...
 '\it{k_{ap}}', '\it{\theta_{ap}}', '\it{\beta_{ap}}'};

% output_all stores full ODE output (Q,P,R) at each timepoint
output_all = zeros(numSamples, length(time_points), 3); % 3: Q, P, R

parfor i = 1:numSamples
    % Create parameter struct for this sample
    params = struct();
    for j = 1:numel(param_names)
        params.(param_names{j}) = param_samples(i, j);
    end

    % Solve ODE with current parameter set
    [~, y] = ode15s(@(t,y) LiverRegenModel_PRCC(t,y,struct2array(params)), time, y0);

    % Initialize temp output for this sample
    temp_output = zeros(length(time_points), 3); % rows = time points, cols = Q,P,R

    % Loop through timepoints and extract Q, P, R
    for k = 1:length(time_points)
        idx = find(time >= time_points(k), 1);
        if isempty(idx)
            warning('No index found for time point %.2f. Skipping...', time_points(k));
            continue;
        end
        if idx > size(y, 1)
            warning('Index %d exceeds the number of rows in y (%d)', idx, size(y, 1));
            continue;
        end
        Q = y(idx, 8);
        P = y(idx, 9);
        R = y(idx, 10);
        temp_output(k, :) = [Q, P, R];
    end

    % Store in final output array
    output_all(i, :, :) = temp_output;
end

% Initialize an array to store PRCC values and p-values for all timepoints
prcc_vals_all = zeros(numel(param_names), length(time_points));  % PRCC values
prcc_pvals_all = zeros(numel(param_names), length(time_points));  % p-values

% Loop through each timepoint for PRCC analysis and plotting
for k = 1:length(time_points)
    curr_timepoint = time_points(k);
    curr_label = timepoint_names{k};

    % Calculate liver volume at this timepoint
    Q = output_all(:, k, 1);
    P = output_all(:, k, 2);
    R = output_all(:, k, 3);
    liver_vol = Q + 2*P + 1.5*R;

    % Rank transform
    X = tiedrank(param_samples);
    LiverVolRank = tiedrank(liver_vol);

    % PRCC computation
    prcc_vals = zeros(numel(param_names), 1);
    prcc_pvals = zeros(numel(param_names), 1);
    for j = 1:numel(param_names)
        other_idx = setdiff(1:numel(param_names), j);
        [~, ~, rx] = regress(X(:,j), [ones(numSamples,1), X(:,other_idx)]);
        [~, ~, ry] = regress(LiverVolRank, [ones(numSamples,1), X(:,other_idx)]);
        [Rmat, Pmat] = corrcoef(rx, ry);
        prcc_vals(j) = Rmat(1,2);
        prcc_pvals(j) = Pmat(1,2);
    end

    % Store PRCC values and p-values for this timepoint
    prcc_vals_all(:, k) = prcc_vals;
    prcc_pvals_all(:, k) = prcc_pvals;
end


