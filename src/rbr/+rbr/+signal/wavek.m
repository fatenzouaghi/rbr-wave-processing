function k = wavek(f, h, g)
% WAVeK  Linear dispersion solver:  w^2 = g k tanh(k h)
% f : frequency [Hz], h : depth [m], g : gravity [m/s^2]
w = 2*pi*f(:);
k = (w.^2./g);               % initial guess (shallow)
for it = 1:30
    kh = k.*h;
    t  = tanh(kh);
    F  = g.*k.*t - w.^2;
    dF = g.*t + g.*k.*(1 - t.^2).*h;
    k  = k - F./max(dF,eps);
end
k = max(k, 0);
end
