function res = compute_block_spectrum(press_samp, depth_mean, opts)
% Compute pressure PSD, convert to elevation PSD (optional attenuation),
% then derive bulk wave parameters 

rho = 1023; g = 9.81;

% 1) Pressure PSD 
[ff, df, PPp, ~, ~] = spectrum (press_samp, opts.fs, opts.nfft, depth_mean);

% 2) Pressure → elevation via linear wave theory
PP = convert_pressure_to_elevation(ff, depth_mean, opts, PPp);

% 3) Frequency window
II = (ff >= opts.minFreq) & (ff <= opts.maxFreq);
ff = ff(II);  PP = PP(II);
if numel(ff) > 1, df = mean(diff(ff)); end

% 4) Spectral moments & bulk parameters
m0 = sum(PP) * df;
m1 = sum(ff .* PP) * df;
m2 = sum((ff.^2) .* PP) * df;

res.Hs   = 4 * sqrt(m0);
[~,pk]   = max(PP);
res.Tp   = 1 ./ ff(pk);
res.Tm01 = m0 / m1;
res.Tm02 = sqrt(m0 / m2);

% 5) IG vs sea–swell partitions
Iig        = ff <  opts.igCutoff;
Isw        = ff >= opts.igCutoff;
res.Hs_IG  = 4 * sqrt(sum(PP(Iig)) * df);
res.Hs_SW  = 4 * sqrt(sum(PP(Isw)) * df);


