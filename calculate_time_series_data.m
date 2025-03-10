function next_timepoint_data = calculate_time_series_data(patient_id, timepoint, initial_values)

% Initial conditions
% All these are from python code
% In python indexing start from 0 so always add 1 for correct index in MATLAB
timepoint = timepoint+1; % in python index starts from 0 and it starts form 1 in MATLAB
patient_id = patient_id;

%% Initial conditions from VAE

[tnf0, jak0, st30, soc0, ecm0, ie0, gf0, q0, p0, r0] = deal(initial_values(1), initial_values(2), initial_values(3), initial_values(4), ...
    initial_values(5), initial_values(6), initial_values(7), initial_values(8), initial_values(9), initial_values(10));

vmp = 2; % volume multiplier for cells in P
vmr = 1.5; % volume multiplier for cells in R

%%  ODE simulation

% Patients' time points(in minutes): Before Surgery | 5 min | 30 min | 60 min | 120 min |
% | 1440 min (1day)| 2880 min (2days)| 4320 min (3days)| 5760 min (4days)|
% | 14400 min (10days)| 129600 min | (3months)| 259200 min (6months)| 525600 min (1year)|

time_points = [1,5,30,60,120,1440,2880,4320,5760,14400,129600,...
    259200,525600]; % time points we require the data | 1 min is taken as Before Suergery for simulation

time = linspace(time_points(timepoint)/60,365*24, (365*24*60 - time_points(timepoint))+1 ); % minute time step

y0 = [tnf0, jak0, st30, soc0, ecm0, ie0, gf0, q0, p0, r0];

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15);
[t,y] = ode15s(@LiverRegenModel, time, y0, options);

% Here if the python timepoint = 7 means matlab timepoint = timepoint + 1 = 8 (Day 3) 
% So, the iteration starts with Day 3 i.e. iteration i = 1 is for time_points(timepoint) = 4320 min (3days)
% ==> So, to get the correct iteration we have to subtract the initial timepoints value here 4320
% Iteration 1 : has the functional value of timepoint 4320 (Day 3)
% Iteration 1440 (we can get by 5760 - 4320): has the functional value of timepoint 5760 min (4days) 
% ===> To get the correct iteration no. --> time_points(timepoint+1:end)-time_points(timepoint)

next_data = transpose(y(time_points(timepoint+1:end)-time_points(timepoint),:)); % timepoint + 1 --> get the next time point data
next_timepoint_data = next_data;

end