function PP = convert_pressure_to_elevation(ff, h_mean, opts, PPp)
% CONVERT_PRESSURE_TO_ELEVATION  Convert pressure spectrum to elevation spectrum
% using linear wave theory with or without attenuation.
%
% Inputs:
%   ff      : frequency vector
%   h_mean  : mean depth approximation
%   opts    : options structure containing parameters
%   PPp     : pressure PSD (Power Spectral Density)
%
% Output:
%   PP      : elevation PSD (Power Spectral Density)

% Physical constants
rho = 1023;  % Density of seawater (kg/m^3)
g = 9.81;    % Gravitational acceleration (m/s^2)

% If attenuation is used, apply the transfer function
if opts.use_attenuation
    % Compute the wave number using linear wave theory
    k = wavek(ff, h_mean, g);
    
    % Compute the transfer function
    TF = cosh(k * h_mean) ./ cosh(k * opts.hd) / (rho * g);  % Transfer function
    
    % Apply the transfer function to the pressure PSD to get elevation PSD
    PP = PPp .* (TF .^ 2);  % Elevation PSD with attenuation
else
    % Without attenuation: simple conversion to elevation PSD
    PP = PPp / (rho * g) ^ 2;
end
end
