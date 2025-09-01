function out = process_rbr_pressure(rskFile, meteoFile, opts)
% PROCESS_RBR_PRESSURE  End-to-end pipeline: RBR pressure → wave parameters.
% - Reads .rsk via RSKtools
% - Barometric correction using meteo .mat (ECCC)
% - Splits the record into fixed-length blocks (delay_sec) and computes spectra
%  via compute_block_spectrum()
% - Pressure→elevation transfer handled in convert_pressure_to_elevation()

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
S = readmatrix(meteoFile); 
Time   = S(:,1);         
Press  = S(:,2);          
Temp     = S(:,3);          

% fill then interpolate onto RBR timestamps
Press_intrp = fillmissing(Press*1000, 'linear'); 
Temp_intrp  = fillmissing(Temp, 'linear');

pAtm = interp1(Time, Press_intrp, time, 'spline');  % Pa
Tmp  = interp1(Time, Temp_intrp,   time, 'spline');  % °C

% Isothermal barometric leveling to sea level
hs         = (R*(Tmp + 273.15)) ./ (M*g);  % scale height (m)
pAtm_corr  = pAtm .* exp(opts.alti ./ hs); % Pa corrected to MSL

% Gauge-like pressure (sensor minus atmosphere)
p = pRaw - pAtm_corr;                       % Pa

%% -------- 3) Block segmentation ----------
ndelay    = max(1, round(opts.delay_sec * opts.fs)); % samples per block
nbspectre = floor(numel(p) / ndelay);

% Allocate outputs
Time        = nan(nbspectre,1);
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
    WL(i)  = pm/(rho*g) + opts.hd;       % m


    % Only compute a spectrum if water depth exceeds threshold
    if WL(i) > opts.crit_m
        res   = compute_block_spectrum(press_samp, WL(i), opts);
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
out.spec.Hs       = Hs;
out.spec.Hs_IG    = Hs_IG;
out.spec.Hs_SW    = Hs_SW;
out.spec.Tp       = Tp;
out.spec.Tm01     = Tm01;
out.spec.Tm02     = Tm02;
end







