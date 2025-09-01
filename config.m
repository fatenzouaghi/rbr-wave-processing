% config.m — default options (no raw data paths, safe to publish)
opts_defaults = struct();

% Sampling & spectra
opts_defaults.fs              = 4;        % Hz
opts_defaults.nfft            = 1024;
opts_defaults.minFreq         = 0.0083;   % Hz
opts_defaults.igCutoff        = 0.05;     % Hz
opts_defaults.maxFreq         = 0.5;      % Hz


% Block length & minimum depth threshold
opts_defaults.delay_sec       = 60*20;    % 20-minute blocks
opts_defaults.crit_m          = 0.35;     % minimum WL to compute spectrum (m)

% Site/sensor
opts_defaults.alti      = 44.5;    % m (station altitude above MSL)
opts_defaults.zmembrane = -1.236;  % m (sensor membrane elevation)
opts_defaults.zbottom   = -1.52;   % m (bed elevation)
