function out = process_rbr_pressure(rskFile, meteoMat, opts)
% PROCESS_RBR_PRESSURE  End-to-end pipeline for RBR pressure → waves.
%
% Inputs
%   rskFile   : path to RBR .rsk
%   meteoMat  : path to meteo .mat with variables:
%               - Time_UTC (datetime, UTC)
%               - Press       (kPa)    [converted to Pa here]
%               - Temperature (°C)
%   opts      : struct with fields (defaults in parentheses):
%               alti (required, m)  – met station altitude above MSL
%               hd  (required, m)  – sensor Height Above Bed = zmembrane - zfond
%               fs   (4)            – sampling frequency [Hz]
%               nfft (1024)         – FFT length
%               use_attenuation (true) – apply pressure→elevation transfer
%               minFreq (0.0083), igCutoff (0.05), maxFreq (0.5)
%
% Output
%   out.Time      : datetime vector (RBR timestamps)
%   out.ABSOLUTE_WL     : absolute water level [m]
 %  out.WL_CGVD2013    : water level [relatif to the CGVD2013, m]
%   out.spec      : struct with spectral results:
%                   Hs, Hs_IG, Hs_SW, Tp, Tm01, Tm02
%
%
% Notes:
% - Reading .rsk requires RSKtools (RBR): see the official page.
% - Meteo data in meteo.mat should come from a reliable source (e.g., ECCC).

% ---- defaults (safe fallbacks; values from config.m override them)
if ~isfield(opts,'fs'),            opts.fs = 4;      end
if ~isfield(opts,'nfft'),          opts.nfft = 1024; end
if ~isfield(opts,'use_attenuation'), opts.use_attenuation = true; end
if ~isfield(opts,'minFreq'),       opts.minFreq = 0.0083;  end
if ~isfield(opts,'igCutoff'),      opts.igCutoff = 0.05;   end
if ~isfield(opts,'maxFreq'),       opts.maxFreq = 0.5;     end


% Physical constants
rho = 1023; g = 9.81; R = 8.314; M = 0.02896;

%% 1) Read RBR time series (absolute pressure) via RSKtools
RSK  = RSKopen(rskFile);
Data = RSKreaddata(RSK);

time = Data.data.tstamp;

% PRESSURE 
pRaw = Data.data.values(:,1) * 1e4;   % dbar → Pa

%% 2) Load meteo.mat (kPa, °C) and interpolate to RBR time
S = load(meteoMat);  % expects Time_UTC, Press (kPa), Temperature (°C)
pAtm = interp1(S.Time_UTC, S.Press*1000, time, 'spline', 'extrap');   % Pa
Temp = interp1(S.Time_UTC, S.Temperature, time, 'spline', 'extrap');  % °C

% Barometric leveling (isothermal approximation)
hs = (R*(Temp + 273.15)) ./ (M*g);         % scale height [m]
pAtm_corr = pAtm .* exp(opts.alti ./ hs);  % Pa adjusted to sea level

% Remove atmospheric pressure → gauge-like pressure
p = pRaw - pAtm_corr;                       % Pa

% Absolute water level: ABSOLUTE_WL = p/(ρg) + hd
ABSOLUTE_WL = p./(rho*g) + opts.hd;

%% 3) Pressure spectrum, then convert to elevation spectrum
%   [ff, df, PP, fp, Hs] = spectrum (e, Fs, nfft, d)
% Pass corrected pressure (Pa), fs, nfft, and a depth proxy d:
h_mean = mean(p)/(rho*g) + opts.hd;   % approximate mean depth used as 'd'
[ff, df, PPp, ~, ~] = spectrum (p, opts.fs, opts.nfft, h_mean); % PPp: pressure PSD

% PRESSURE → ELEVATION (linear wave theory)
if opts.use_attenuation
    k  = wavek(ff, h_mean, g);
    TF = cosh(k.*h_mean) ./ cosh(k.*opts.hd) / (rho*g)
    PP = PPp .* (TF.^2);                        % elevation PSD [m^2/Hz]
else
    PP = PPp ./ (rho*g).^2;                     % hydrostatic-only fallback
end

% Frequency band selection
II = (ff >= opts.minFreq & ff <= opts.maxFreq);
ff = ff(II);  PP = PP(II);

% Robust df (if frequency grid not strictly uniform)
if numel(ff) > 1
    df = mean(diff(ff));
end

%% 4) Spectral moments → bulk wave parameters
m0 = sum(PP)          * df;
m1 = sum(ff .* PP)    * df;
m2 = sum((ff.^2) .* PP) * df;

Hs   = 4*sqrt(m0);
[~,pk] = max(PP);  Tp = 1./ff(pk);
Tm01 = m0 / m1;
Tm02 = sqrt(m0 / m2);

% IG / sea–swell split
Iig  = ff <  opts.igCutoff;
Isw  = ff >= opts.igCutoff;
HsIG = 4*sqrt(sum(PP(Iig)) * df);
HsSW = 4*sqrt(sum(PP(Isw)) * df);

%% 5) Pack outputs
out.Time        = Time(:);
out.ABSOLUTE_WL       = ABSOLUTE_WL(:);
out.WL_CGVD2013      = WL_CGVD2013 (:);
out.spec.Hs     = Hs;
out.spec.Hs_IG   = Hs_IG;
out.spec.Hs_SW   = Hs_SW;
out.spec.Tp     = Tp;
out.spec.Tm01   = Tm01;
out.spec.Tm02   = Tm02;

end





