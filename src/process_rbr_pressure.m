function out = process_rbr_pressure(rskFile, meteoMat, opts)
% PROCESS_RBR_PRESSURE  End-to-end pipeline for RBR pressure → waves.
%
% This function processes the raw pressure data from an RBR sensor and 
% calculates wave parameters (Hs, Tp, Tm01, Tm02) from the pressure spectra 
% and meteorological data. 
%
% Inputs
%   rskFile   : path to the .rsk file (RBR sensor data)
%   meteoMat  : path to the .mat file containing meteorological data with the 
%               following variables:
%               - Time_UTC (datetime, UTC)
%               - Press (Pressure in kPa)
%               - Temperature (Temperature in °C)
%   opts      : structure with parameters, including the following fields 
%               (defaults are listed in parentheses):
%               alti (m)              : meteorological station altitude above sea level (required)
%               zmembrane (m)         : Sensor height above the bed (zmembrane) (required)
%               zbottom (m)           : Depth at the bottom (zbottom) (required)
%               fs   (4)               : Sampling frequency in Hz [optional]
%               nfft (1024)            : FFT size [optional]
%               use_attenuation (true): Apply pressure → elevation transfer [optional]
%               minFreq (0.0083)      : Minimum frequency for spectrum [optional]
%               igCutoff (0.05)       : Infragravity wave frequency cutoff [optional]
%               maxFreq (0.5)         : Maximum frequency for spectrum [optional]
%
% Outputs
%   out.Time       : datetime vector of RBR timestamps
%   out.ABSOLUTE_WL: Absolute water level in meters
%   out.WL_CGVD2013: Water level relative to CGVD2013 in meters
%   out.spec       : structure with spectral results:
%                    Hs, Hs_IG, Hs_SW, Tp, Tm01, Tm02
%
% Notes:
% - RSKtools is required to read .rsk files (see official RBR page).
% - The meteo.mat file should contain reliable meteorological data (e.g., from ECCC).

% ---- Load the configuration file with parameters
run('path_to_config/config.m');  % Adjust path to where your config.m file is located

% If opts is provided, merge with default options from config.m
if nargin > 2
    opts = merge_opts(config_opts, opts);  % Merge any provided options with defaults
end

% Calculate hd (sensor height above the bed)
opts.hd = opts.zmembrane - opts.zbottom;  % Sensor height above the bed (zmembrane - zbottom)

% Physical constants
rho = 1023;  % Density of seawater (kg/m^3)
g = 9.81;    % Gravitational acceleration (m/s^2)

%% Set Additional Parameters (delay and crit)
delay = 60 * 20;  % Time interval for spectrum calculation (20 minutes)
crit = 0.35;  % Minimum height for which the spectrum calculation will be done

%% 1) Read RBR data (.rsk) using RSKtools
% Load the RBR sensor data using RSKtools
RSK = RSKopen(rskFile);
Data = RSKreaddata(RSK);
time = Data.data.tstamp;

% Raw pressure data (dbar to Pa conversion)
pRaw = Data.data.values(:,1) * 1e4;   % Convert from dbar to Pa

%% 2) Load and interpolate meteorological data to RBR timestamps
% The .mat file must contain variables: Time_UTC, Press (in kPa), Temperature (°C)
S = load(meteoMat);  % Load meteorological data

% Interpolate meteorological data to match RBR timestamps
pAtm = interp1(S.Time_UTC, S.Press * 1000, time, 'spline', 'extrap');  % Convert pressure to Pa
Temp = interp1(S.Time_UTC, S.Temperature, time, 'spline', 'extrap');   % Temperature in °C

% Barometric leveling using isothermal approximation
hs = (R * (Temp + 273.15)) ./ (M * g);         % Scale height in meters
pAtm_corr = pAtm .* exp(opts.alti ./ hs);       % Atmospheric pressure corrected to sea level

% Subtract atmospheric pressure to get gauge pressure (relative to sea level)
p = pRaw - pAtm_corr;                         % Pressure in Pa

% Calculate absolute water level (in meters): H = p / (ρ * g) + hd
ABSOLUTE_WL = p / (rho * g) + opts.hd;


%% 4) Spectral Calculation Loop (Using delay and crit)

% Calculate the number of spectra based on the data length and delay
nbspectre = floor(length(p) / delay); % Number of spectra

Hs = nan(nbspectre, 1);  % Initialize significant wave height array
Time = nan(nbspectre, 1); % Initialize time array
ABSOLUTE_WL = Time;  % Initialize water level array

% Loop through each spectrum
for i = 1:nbspectre
    disp(['Processing spectrum ', num2str(i), ' / ', num2str(nbspectre)])

    % Select data for the current spectrum (based on delay)
    press_samp = p((i - 1) * delay + 1:i * delay);   
    Time(i) = mean(time((i - 1) * delay + 1:i * delay));  % Average time of the spectrum

    pm = mean(press_samp);
    ABSOLUTE_WL(i) = pm / (rho * g);  % Absolute water level in Pascal
    ABSOLUTE_WL(i) = ABSOLUTE_WL(i) + opts.hd;  % Adjusted absolute water level

    % **Apply the threshold check (crit)**
    if ABSOLUTE_WL(i) > crit
        % Perform spectrum calculation (only once)
        [ff, df, PP, ~, ~] = spectrum (press_samp, opts.fs, opts.nfft, ABSOLUTE_WL(i));
        km = wavek(ff, ABSOLUTE_WL(i), g);  % Wave number

        % Apply attenuation and convert to elevation
        PP = convert_pressure_to_elevation(ff, ABSOLUTE_WL(i), opts, PP); 
        
        % Calculate Hs, Tp, Tm01, etc., here
        Hs(i) = 4 * sqrt(sum(PP * df));  % Significant wave height
        % Continue with other parameters (Tp, Tm01, Tm02) here...
    end
end

%% 5) Pack results into the output structure
out.Time = Time(:);
out.ABSOLUTE_WL = ABSOLUTE_WL(:);
out.WL_CGVD2013 = ABSOLUTE_WL + opts.zbottom;  % Water level relative to CGVD2013
out.spec.Hs = Hs;
% Continue packing other parameters (Hs_IG, Hs_SW, Tp, Tm01, Tm02) here...

end





