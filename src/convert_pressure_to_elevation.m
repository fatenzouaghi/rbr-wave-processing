function PP = convert_pressure_to_elevation(ff, depth_mean, opts, PPp)

rho = 1023; 
g   = 9.81;

k  = wavek(ff, depth_mean, g);                         % wave number
TF = cosh(k.*depth_mean) ./ cosh(k.*opts.hd) / (rho*g);% pressure→elevation 
PP = PPp .* (TF.^2);                                   % elevation PSD [m^2/Hz]
end
