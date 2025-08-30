function out = process_rbr_pressure(rskFile, meteoMat, opts)
% PROCESS_RBR_PRESSURE  End-to-end pipeline: RBR pressure → wave parameters.
% - Reads .rsk via RSKtools
% - Barometric correction using meteo .mat (ECCC)
% - Splits the record into fixed-length blocks (delay_sec) and computes spectra
%   per block via compute_block_spectrum()
% - Pressure→elevation transfer handled in convert_pressure_to_elevation()
%
% Inputs
%   rskFile   : path to RBR .rsk file
%   meteoCSV  : path to .csv with variables:
%               - Time_UTC (datetime, UTC)
%               - Press (kPa)
%               - Temperature (°C)
%   opts      : struct overriding defaults from config/config.m
%
% Output (struct)
%   out.Time, out.ABSOLUTE_WL, out.WL_CGVD2013
%   out.spec = struct('Hs','Hs_IG','Hs_SW','Tp','Tm01','Tm02')

%% -------- 0) Load defaults & merge user opts ----------
run('config.m');                 % defines opts_defaults
if nargin < 3 || isempty(opts)
    opts = opts_defaults;
else
    opts = merge_opts(opts_defaults, opts);
end

% Sensor-derived parameter
opts.hd = opts.zmembrane - opts.zbottom; % sensor height above bed (m)

% Physical constants
rho = 1023;               % seawater density (kg/m^3)
g   = 9.81;               % gravity (m/s^2)
R   = 8.314;              % gas constant (J/(mol·K))
M   = 0.02896;            % molar mass of air (kg/mol)

%% -------- 1) Read RBR (.rsk) ----------
RSK  = RSKopen(rskFile);
Data = RSKreaddata(RSK);
time = Data.data.tstamp;            % datetime (UTC)
pRaw = Data.data.values(:,1) * 1e4; % dbar → Pa (absolute pressure)

%% -------- 2) Load meteo & barometric leveling ----------
S = load(meteoCSV);  % expects: Time_UTC, Press (kPa), Temperature (°C)
Time_UTC   = S.Time_UTC;         % datetime (UTC)
Press_kPa  = S.Press;            % kPa
Temp_C     = S.Temperature;      % °C

% Safe fill then interpolate onto RBR timestamps
Press_Pa = fillmissing(S.Press*1000, 'linear'); % kPa→Pa
Temp_C   = fillmissing(S.Temperature, 'linear');

pAtm = interp1(S.Time_UTC, Press_Pa, time, 'spline', 'extrap');  % Pa
Tmp  = interp1(S.Time_UTC, Temp_C,   time, 'spline', 'extrap');  % °C

% Isothermal barometric leveling to sea level
hs         = (R*(Tmp + 273.15)) ./ (M*g);  % scale height (m)
pAtm_corr  = pAtm .* exp(opts.alti ./ hs); % Pa corrected to MSL

% Gauge-like pressure (sensor minus atmosphere)
p = pRaw - pAtm_corr;                       % Pa

%% -------- 3) Block segmentation ----------
ndelay    = max(1, round(opts.delay_sec * opts.fs)); % samples per block
nbspectre = floor(numel(p) / ndelay);

% Allocate outputs
Time        = NaT(nbspectre,1);
ABSOLUTE_WL = nan(nbspectre,1);
WL_CGVD2013 = nan(nbspectre,1);
Hs          = nan(nbspectre,1);
Hs_IG       = nan(nbspectre,1);
Hs_SW       = nan(nbspectre,1);
Tp          = nan(nbspectre,1);
Tm01        = nan(nbspectre,1);
Tm02        = nan(nbspectre,1);

%% -------- 4) Loop blocks & call spectral function ----------
for i = 1:nbspectre
    i1 = (i-1)*ndelay + 1;
    i2 = i*ndelay;

    press_samp = p(i1:i2);
    Time(i)    = mean(time(i1:i2));

    % Absolute water level at the sensor for this block (m)
    pm              = mean(press_samp);           % Pa
    ABSOLUTE_WL(i)  = pm/(rho*g) + opts.hd;       % m
    WL_CGVD2013(i)  = ABSOLUTE_WL(i) + opts.zbottom;

    % Only compute a spectrum if water depth exceeds threshold
    if ABSOLUTE_WL(i) > opts.crit_m
        res   = compute_block_spectrum(press_samp, ABSOLUTE_WL(i), opts);
        Hs(i)    = res.Hs;
        Hs_IG(i) = res.Hs_IG;
        Hs_SW(i) = res.Hs_SW;
        Tp(i)    = res.Tp;
        Tm01(i)  = res.Tm01;
        Tm02(i)  = res.Tm02;
    end
end

%% -------- 5) Package output ----------
out.Time          = Time(:);
out.ABSOLUTE_WL   = ABSOLUTE_WL(:);
out.WL_CGVD2013   = WL_CGVD2013(:);
out.spec.Hs       = Hs;
out.spec.Hs_IG    = Hs_IG;
out.spec.Hs_SW    = Hs_SW;
out.spec.Tp       = Tp;
out.spec.Tm01     = Tm01;
out.spec.Tm02     = Tm02;
end





