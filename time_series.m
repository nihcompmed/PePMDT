%% Liver Regeneration Simulation
% This is a time series simulation of the model system with the patient
% data: we have used the initial liver mass of the patients in the model
% and simulate the time series of the 10 model variables

% The data then is used in to neural networks to find the gene clusters
% that represent the model variables

%% Patients data
load liver_resected.mat % percentage of liver resected for 12 patients

% We need the data in the fraction and q0 is the initial fraction remaining

q_zero = 1-(liver_resected/100); % initial liver fractions for 12 patients

%% Initial conditions

tnf0 = 1;
jak0 = 1;
st30 = 1;
soc0 = 1;
ecm0 = 1;
ie0 = 1;
gf0 = 1;
% select q0 based on the patient
% q0 = 0.05; % ~0.35 fraction based on 12 patients data
p0 = 0;
r0 = 0;

vmp = 2; % volume multiplier for cells in P
vmr = 1.5; % volume multiplier for cells in R

%%  ODE simulation
% simulate the model to generate patient wise data

% time = linspace(0,365*24,8761); % hour time step
time = linspace(0,365*24,525601); % minute time step
% Patients' time points(in minutes): | 5 min | 30 min | 60 min | 120 min |
% | 180 min (3hrs)| 240 min (4hrs)| 1440 min (1day)| 2880 min (2days)|
% | 4320 min (3days)| 5760 min (4days)| 14400 min (10days)| 129600 min
% (3months)| 259200 min (6months)| 525600 min (1year)|

% time_points = [5,30,60,120,180,240,1440,2880,4320,5760,14400,129600,...
%     259200,525600]; % time points we require the data

time_points = [1,5,30,60,120,1440,2880,4320,5760,14400,129600,...
    259200,525600]; % time points we require the data
% removed 3hrs and 4hrs as these time points are missing from every patient
% these time points are only mentioned in the article


for i = 1:length(q_zero)
    q0 = q_zero(i); % initial liver fraction patient wise
    y0 = [tnf0, jak0, st30, soc0, ecm0, ie0, gf0, q0, p0, r0];
    [t,y] = ode15s(@LiverRegenModel, time, y0);

    % Get the patient wise time points and save model simulated data only for
    % those time points
    eval(['simulated_data_p',num2str(i),' = transpose(y(time_points+1,:));']) % first iteration is initial values so take from the 2nd
    eval(['save simulated_data_p',num2str(i),' simulated_data_p',num2str(i)])
end

%% Figure: Patient wise time series of model state variables

noRows = 3;
noColumns = 4;

fig1 = figure;
fig1.WindowState = 'maximized';
for i = 1:numel(q_zero)
    subplot(noRows,noColumns,i)
    eval(['data = simulated_data_p',num2str(i),';'])
    plot(1:length(time_points),data','-o','markersize',8,'LineWidth',2.5)

    ax = gca;
    ax.FontSize = 16;
    ax.LineWidth = 1.8;

    cmap = colormap(parula(10)); % colormap for 10 variables
    ax.ColorOrder = cmap;
    % mylinestyles = ["-"; "--"; "-o"];
    % ax.LineStyleOrder = mylinestyles;
    
    ax.XLim = [0.5 length(time_points)+0.5];
    ax.XTick= 1:1:length(time_points);
    % ax.XTickLabel = {'5mi','30mi','60mi','120mi','3hrs','4hrs',...
    %     '1day','2days','3days','4days','10days','3mon','6mon','1year'};
    ax.XTickLabel = {'1min','5mi','30mi','60mi','120mi','1day','2days','3days',...
        '4days','10days','3mon','6mon','1year'};

    title(strcat('Patient',{' '},num2str(i)));

    % set(gca,'FontWeight','bold', 'FontName','Calibri')
    set(gca,'TickLength',[0.005, 0.01])
end

variables = {'TNF','JAK','STAT3','SOCS3','ECM','IE','GF','Q','P','R'};
legend(variables) % show legend for only one plot

% Indicate the high value curves
text(8,data(2,6),"JAK",'Color',cmap(2,:),'FontSize',12);
text(8,data(4,6),"SOCS3",'Color',cmap(4,:),'FontSize',12);


