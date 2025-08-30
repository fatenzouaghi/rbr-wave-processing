function PP = convert_pressure_to_elevation(ff, depth_mean, opts, PPp)
% CONVERT_PRESSURE_TO_ELEVATION
% Convert pressure spectrum (PPp) to surface elevation spectrum (PP)

% where h = depth_mean, z_sensor = opts.hd (height above bed).
rho = 1023; g = 9.81;

if opts.use_attenuation
    k  = wavek(ff, depth_mean, g);                   % user-provided function
    TF = cosh(k.*depth_mean) ./ cosh(k.*opts.hd) / (rho*g);
    PP = PPp .* (TF.^2);
else
    PP = PPp ./ (rho*g).^2;
end
end

