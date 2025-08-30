
% Minimal reproducible example (no figures).
% Requires: MATLAB + RSKtools (RBR) installed and on the path.

% --- Locate repo root (this file sits in examples/) and add src/ to MATLAB path
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));   % add src/ and all its subfolders

% --- Site parameters for THIS example
HAB  = 0.406;   % m  (sensor Height Above Bed = zmembrane - zbottom)
alti = 44.5;    % m  (meteorological station altitude above MSL)

% --- Processing options 
opts = struct( ...
  'alti',alti,'HAB',HAB, ...
  'fs',4,'nfft',1024, ...
  'use_attenuation',true, ...
  'minFreq',0.0083,'igCutoff',0.05,'maxFreq',0.5,'transferStop',0.40);

% --- Example files provided in data/
rskFile  = fullfile(root,'data','206599_20220715_2119.rsk');
meteoMat = fullfile(root,'data','meteo.mat');   % variables: Time_UTC (datetime, UTC), Press (kPa), Temperature (°C)

% --- Run the end-to-end processing
out = process_rbr_pressure(rskFile, meteoMat, opts);

% --- Numeric summary (includes water level stats)
L = out.level; L = L(~isnan(L));  % ignore NaNs for stats
summaryRow = table(out.Hs, out.HsIG, out.HsSW, out.Tp, out.Tm01, out.Tm02, ...
                   mean(L,'omitnan'), min(L), max(L), ...
                   'VariableNames', {'Hs','HsIG','HsSW','Tp','Tm01','Tm02','MeanLevel','MinLevel','MaxLevel'});
disp(summaryRow)

% --- (Optional) Export time–level series as CSV (no plotting)
TT = timetable(out.time, out.level, 'VariableNames', {'level'});
writetimetable(TT, fullfile(fileparts(mfilename('fullpath')), 'level_example.csv'));

% End of example
