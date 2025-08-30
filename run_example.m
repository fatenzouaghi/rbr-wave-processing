% Add paths and run the pipeline on one deployment
addpath('src');  % code
run('config/config.m');  % defines opts_defaults

% ---- user inputs (raw data live on Zenodo) 
rskFile  = 'data/raw/rbr_file.rsk';
meteoMat = 'data/meteo/meteo_file.mat'; % contains Time_UTC, Press(kPa), Temperature(C)

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
