
%% PART 1 - DATA IMPORT & CLEANING

% import data
importdata;
figure_casttimes;
figure_replicates;

% adjust sample times so that they are uniform
uniformsampling;
figure_uniformsampling;
figure_uniformsampling_error;

% split into separate datasets by depth and omit contigs that are absent
splitdepths;
interpolatenans; % interpolate for missing values (use caution!)
figure_timeseries;
figure_timeseries_summary;

fprintf('PART 1 COMPLETE\n\n');


%% PART 2 - PARTIAL AUTOCORRELATION & AUTOREGRESSION

% partial autocorrelation
partialautocorrelation;
figure_partialautocorrelation;
figure_partialautocorrelation_summary;

% autoregression
autoregression;
figure_autoregression;
figure_autoregression_summary;

fprintf('PART 2 COMPLETE\n\n');


%% PART 3 - REGRESSION

% regression
regression;
figure_regression;
figure_regression_summary;

fprintf('PART 3 COMPLETE\n\n');
