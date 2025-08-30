
% config.m - Configuration file for wave analysis parameters

% Sampling frequency (Hz)
opts.fs = 4;

% FFT size
opts.nfft = 1024;

% Use attenuation for pressure to elevation conversion (if applicable) 
opts.use_attenuation = true;

% Frequency bounds
opts.minFreq = 0.0083;   % Minimum frequency
opts.igCutoff = 0.05;    % Infragravity wave cutoff frequency
opts.maxFreq = 0.5;      % Maximum frequency

% Sensor parameters
opts.alti = 44.5;   % Altitude of the meteorological station above sea level (meters)
opts.zmembrane = -1.236;   % Height of the sensor above the bed (meters)
opts.zbottom = -1.52; % Depth at the bottom (meters)
