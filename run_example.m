% Add paths and run the pipeline on one deployment
addpath('src');  % code
run('config/config.m');  % defines opts_defaults

% ---- user inputs (raw data live on Zenodo) 
% --- EDIT HERE: put your LOCAL paths (examples below) ---------------------
% Windows example:
% rskFile  = 'C:\Users\me\Downloads\1_1_206715_20210906_2304.rsk';
% meteoMat = 'C:\Users\me\Downloads\METEO_GF_2021_UTC.mat';

% macOS example:
% rskFile  = '/Users/me/Downloads/1_1_206715_20210906_2304.rsk';
% meteoMat = '/Users/me/Downloads/METEO_GF_2021_UTC.mat';

% Linux example:
% rskFile  = '/home/me/Downloads/1_1_206715_20210906_2304.rsk';
% meteoMat = '/home/me/Downloads/METEO_GF_2021_UTC.mat';

% >>> Replace the two lines below with YOUR local paths:
rskFile  = '/ABSOLUTE/OR/RELATIVE/PATH/TO/your_file.rsk';
meteoMat = '/ABSOLUTE/OR/RELATIVE/PATH/TO/METEO_GF_2021_UTC.mat';

% Optional: override a few defaults at runtime
opts = struct();
% opts.fs = 4; % keep default unless you need to override

% ---- call main processing ----
out = process_rbr_pressure(rskFile, meteoMat, opts);

% ---- export CSV (drop NaN/NaT rows) ----
if ~exist('output','dir'), mkdir('output'); end
csvPath = fullfile('output', 'waves_timeseries.csv');
export_csv(out, csvPath);
fprintf('Wrote: %s\n', csvPath);
