% Script to select a representative period from a timeseries of electricity
% prices in Austria

scenario = 'Leviathan';
year = '2035';

ElData_Weekly = readmatrix(strcat('Timeseries_data/El_price_',scenario,'/',scenario,'_',year,'_weekly.csv'));
ElData_Annual = readmatrix(strcat('Timeseries_data/El_price_',scenario,'/',scenario,'_',year,'.csv'));

%% Parameters for representative period selection
n_total = 13; % number of periods the full data is partitioned to
n_repr = 1; % number of representative periods to choose from full data
n_bin = 50; % number of bins to partition the DCs for optimization
max_value = 530;
% Max Shinyhappy 2035: 350
% Max Shinyhappy 2045: 500
% Max Leviathan 2035: 530
% Max Leviathan 2045: 670

% Dividing duration curves into bins
bins = linspace(0,max_value-(max_value/n_bin),n_bin);
for i = 0:3
    % Parameter Lpb: The share of time during which the original time series p
    % has a value greater than or equal to the lowest value in the range
    % corresponding to bin b
    AnnualData = ElData_Annual(i*2184+1:i*2184+2184,:);
    DataEl = ElData_Weekly(:,i*13+1:i*13+13);
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
    
    Apbd_acc = zeros(length(DataEl),n_bin,n_total);
    for h = 1:n_total
    Apbd_acc(:,:,h) = bsxfun(@ge, repmat(DataEl(:,h),1,n_bin),repmat(bins,length(DataEl),1));
    end
    
    % Convert the logical matrix to double (1 for true, 0 for false)
    Apbd_acc = double(Apbd_acc);
    Apbd = sum(Apbd_acc)/length(DataEl);
    Apbd = reshape(Apbd,n_bin,n_total)';
    
    %% Optimization
    % Initialize Constraints and Objective function 
    Con = [];
    obj = 0;
    
    % Variables
    error = sdpvar(1,n_bin,'full');                                  % variable for error
    wd = sdpvar(n_total,1,'full');                                   % weights per potential representative time period
    ud = binvar(n_total,1,'full');                                   % binary variable: repr. period selected
    
    
    % constraints
    Con = [Con, error == abs(Lpb-sum(repmat((wd/n_total),1,n_bin).*Apbd)),            % Error function
    sum(ud) == n_repr,                                           % Number of selected periods
    sum(wd) == n_total,                                          % sum of weights is n_total
    wd <= ud * n_total,                                          % Every weight in wd has to be smaller than n_total, if the period is selected
    wd >= 0,
    ];
    
    
    %% Objective function
    obj = obj + sum(error,"all");
    
    %% optimize
    ops = sdpsettings('verbose',1,'debug',1);
    optimize(Con,obj,ops)
    
    %weights = zeros(13,4)
    %selectedPer = zeros(13,4)
    weights(:,i+1) = value(wd);
    selectedPer(:,i+1) = value(ud);
    Periods(1+n_total*i:n_total+n_total*i,:) = value(ud);
end
%% Post Optimization
Weeks = Periods';
Weeksplusempty = ElData_Weekly.*repmat(Weeks,length(ElData_Weekly),1);
Weeksplusempty(isnan(Weeksplusempty)) = 0;
% Find columns with all zeros
zeroColumns = all(Weeksplusempty == 0, 1);
    
% Remove zero columns
SelectedWeeks= Weeksplusempty(:, ~zeroColumns);
Weekscolumn = SelectedWeeks(:);
writematrix(Weekscolumn,strcat('Timeseries_data/El_price_',scenario,'/',scenario,'_',year,'_4Periods.csv'));