function out = process_rbr_pressure(rskFile, meteoMat, opts)
% PROCESS_RBR_PRESSURE  End-to-end RBR processing (pressure → level + spectra)

arguments
  rskFile (1,1) string
  meteoMat (1,1) string
  opts.alti (1,1) double
  opts.HAB  (1,1) double
  opts.fs   (1,1) double = 4
  opts.use_attenuation (1,1) logical = true
  opts.minFreq double = 0.0083
  opts.igCutoff double = 0.05
  opts.maxFreq double = 0.5
  opts.transferStop double = 0.40
end

% --- read rsk and meteorlogical timeseries
[time,pRaw] = read_rsk_timeseries(rskFile);           % Time (UTC), Pressure (dbar)
pRaw = pRaw * 1e4;   % dbar → Pa
[S_t,S_p_kPa,S_T] = load_meteo_mat(meteoMat);         % Time (UTC), atmospheric pressure (kPa) , temperature ( °C)

% --- meteo interp + kPa→Pa
pAtm = interp1(S_t, S_p_kPa*1000, time, 'spline', 'extrap');
Temp = interp1(S_t, S_T,        time, 'spline', 'extrap');

% --- barometric leveling
g=9.81; R=8.314; M=0.02896;
hs = (R*(Temp+273.15))./(M*g);
pAtm_corr = pAtm .* exp(opts.alti ./ hs);

% --- remove Patm and compute level
rho=1023;
p   = pRaw - pAtm_corr;               % Pa
level = p./(rho*g) + opts.HAB;        % m

% --- spectra 
spec = spectrum_fabrice(p, time, struct( ...
  'fs', opts.fs, 'HAB', opts.HAB, 'use_attenuation', opts.use_attenuation, ...
  'minFreq', opts.minFreq, 'igCutoff', opts.igCutoff, ...
  'maxFreq', opts.maxFreq, 'transferStop', opts.transferStop));

% --- output
out.time  = time(:);
out.level = level(:);
out.spec  = spec;
end
