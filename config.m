function cfg = config()
% CONFIG  Central configuration for rbr-wave-processing.
%
% References:
% - RSKtools (RBR) to read .rsk files: https://rbr-global.com/support/matlab-tools/
% - Environment and Climate Change Canada (ECCC) – historical met data:
%   https://climat.meteo.gc.ca/historical_data/search_historic_data_f.html
%
% Usage:
%   addpath(genpath("src"));
%   cfg  = config();
%   opts = struct('alti',cfg.alti,'HAB',cfg.HAB,'fs',cfg.fs,...
%                 'use_attenuation',cfg.use_attenuation);
%   out  = process_rbr_pressure(cfg.rskFile, cfg.meteoMat, opts);

%% 1) FILES (edit these paths if needed)
cfg.rskFile  = fullfile("data","206599_20220715_2119.rsk");  % example RBR .rsk
cfg.meteoMat = fullfile("data","meteo.mat");                 % expects Time_UTC, Press(kPa), Temperature(°C)

%% 2) SITE / GEOMETRY (same vertical datum for zmembrane & zfond)
cfg.zmembrane = -1.191;     % m
cfg.zfond     = -1.597;     % m
cfg.HAB       = cfg.zmembrane - cfg.zfond;  % m (Height Above Bed)
cfg.alti      = 44.5;       % m (met station altitude above MSL) — change for your station

%% 3) ACQUISITION / PROCESSING
cfg.fs    = 4;              % Hz sampling frequency
cfg.nfft  = 1024;           % FFT length (Δf = fs/nfft ≈ 0.0039 Hz)
cfg.delay = 60*20;          % s, step between consecutive spectra (20 min)
cfg.crit  = 0.2;           % m, “sensor under water” threshold (if used)

%% 4) FREQUENCY BANDS
cfg.minFreq      = 0.0083;   % Hz 
cfg.igCutoff     = 0.05;    % Hz (IG / sea-swell split)
cfg.maxFreq      = 0.5;     % Hz (upper useful band)
cfg.transferStop = 0.40;    % Hz (cap pressure→elevation transfer above this)

%% 5) PHYSICAL CONSTANTS
cfg.rho = 1023;             % kg/m^3 (seawater)
cfg.g   = 9.81;             % m/s^2
cfg.R   = 8.314;            % gas constant
cfg.M   = 0.02896;          % molar mass of air (kg/mol)

%% 6) OPTIONS
cfg.use_attenuation = true; % apply pressure→elevation transfer function
end
