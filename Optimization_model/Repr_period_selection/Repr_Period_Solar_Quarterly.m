% Script to select a representative period from a timeseries of solar
% irradiation at the location in Judenburg, Austria

%ElData = readmatrix('Solar timeseries/Solar_Weekly.csv');
%QuarterlyData = readmatrix('Solar timeseries/Timeseries_1kWp_2020.csv');
fileWeekly = "Solar timeseries/Solar_quarterly/Solar_Weekly_Quarterly_2020.xlsx";
fileQuarterly = "Solar timeseries/Solar_quarterly/Solar_Quarterly_2020.xlsx";

Quarter = 4;

sheetNames = ['Q',num2str(Quarter)];
Data= readmatrix(fileWeekly,'Sheet',sheetNames);
AnnualData_all = readmatrix(fileQuarterly);
AnnualData = AnnualData_all(:,Quarter);
% calculate duration curves for all time series
n_series = 1;
for i = 1:n_series
    DC(:,i) = sort(AnnualData(:,i),'descend'); 
end

%% Parameters for representative period selection
n_total = [13,13,13,12]; % number of periods the full data is partitioned to
n_repr = 1; % number of representative periods to choose from full data
n_bin = 20; % number of bins to partition the DCs for optimization

% Dividing duration curves into bins
bins = linspace(0,950,n_bin);

% Parameter Lpb: The share of time during which the original time series p
% has a value greater than or equal to the lowest value in the range
% corresponding to bin b

Lpb_acc = zeros(length(AnnualData),n_bin);
% Condition, checking whether the value in AnnualData fits into bin
conditionMatrix = bsxfun(@ge, AnnualData, repmat(bins,length(AnnualData),1));

% Convert the logical matrix to double (1 for true, 0 for false)
Lpb_acc = double(conditionMatrix);
Lpb = sum(Lpb_acc, 1) / length(AnnualData);

% Calculate Apbd
% Parameter Apbd: The share of time during which a potential representative
% time period d of the time series p has a value greater than or equal to 
% the lowest value in the range corresponding to bin b

Apbd_acc = zeros(length(Data),n_bin,n_total(Quarter));
for h = 1:n_total(Quarter)
    Apbd_acc(:,:,h) = bsxfun(@ge, repmat(Data(:,h),1,n_bin),repmat(bins,length(Data),1));
end

% Convert the logical matrix to double (1 for true, 0 for false)
Apbd_acc = double(Apbd_acc);
Apbd = sum(Apbd_acc)/length(Data);
Apbd = reshape(Apbd,n_bin,n_total(Quarter))';

%% Optimization
% Initialize Constraints and Objective function 
Con = [];
obj = 0;

% Variables
error = sdpvar(1,n_bin,'full');                                  % variable for error
wd = sdpvar(n_total(Quarter),1,'full');                                   % weights per potential representative time period
ud = binvar(n_total(Quarter),1,'full');                                   % binary variable: repr. period selected


% constraints
Con = [Con, error == abs(Lpb-sum(repmat((wd/n_total(Quarter)),1,n_bin).*Apbd)),            % Error function
    sum(ud) == n_repr,                                           % Number of selected periods
    sum(wd) == n_total(Quarter),                                          % sum of weights is n_total
    wd <= ud * n_total(Quarter),                                          % Every weight in wd has to be smaller than n_total, if the period is selected
    wd >= 0,
];


%% Objective function
obj = obj + sum(error,"all");

%% optimize
ops = sdpsettings('verbose',1,'debug',1);
optimize(Con,obj,ops)

%weights = zeros(13,4)
%selectedPer = zeros(13,4)
weights(1:n_total(Quarter),Quarter) = value(wd)
selectedPer(1:n_total(Quarter),Quarter) = value(ud)

