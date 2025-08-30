% Add paths and run the pipeline on one deployment
addpath(genpath('src'));     % ajoute src/ et ses sous-dossiers au path
which process_rbr_pressure -all
run('config.m');  % defines opts_defaults

% ---- user inputs (raw data live on Zenodo) 
% --- EDIT HERE: put your LOCAL paths (examples below) ---------------------
% Windows example:
% rskFile  = 'C:\Users\me\Downloads\1_1_206715_20210906_2304.rsk';
% meteoCSV = 'C:\Users\me\Downloads\METEO_GF_2021_UTC.csv';

% macOS example:
% rskFile  = '/Users/me/Downloads/1_1_206715_20210906_2304.rsk';
% meteoCSV = '/Users/me/Downloads/METEO_GF_2021_UTC.csv';

% Linux example:
% rskFile  = '/home/me/Downloads/1_1_206715_20210906_2304.rsk';
% meteoFile = '/home/me/Downloads/METEO_GF_2021_UTC.csv';

% >>> Replace the two lines below with YOUR local paths:
rskFile  = 'Y:\25_NU_COAST\02_GF\01_DATA\02_RBR\2021\01_RAW_DATA\P1.1\1_1_206715_20210906_2304.rsk';
meteoFile = 'Y:\25_NU_COAST\02_GF\01_DATA\02_RBR\2021\03_ATMO\METEO_GF_2021_UTC.csv';

% Optional: override a few defaults at runtime
opts = struct();
% opts.fs = 4; % keep default unless you need to override

% ---- call main processing ----
out = process_rbr_pressure(rskFile, meteoFile, opts);

% ---- export CSV 
T = table( ...
    out.Time, ...
    out.spec.Hs, out.spec.Hs_IG, out.spec.Hs_SW, ...
    out.spec.Tp, out.spec.Tm01, out.spec.Tm02, ...
    out.ABSOLUTE_WL, out.WL_CGVD2013, ...
    'VariableNames', {'Time','HS','Hs_IG','Hs_SW','Tp','Tm01','Tm02','ABSOLUTE_WL','WL_CGVD2013'} );


if ~exist('output','dir'), mkdir('output'); end
writetable(T, fullfile('output','waves_timeseries.csv'));
