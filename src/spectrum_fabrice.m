function [ff, df, PP] = spectrum (e, Fs, nfft)
% Reproduces your windowing/overlap/Hanning approach in a compact way.
% Outputs one-sided PSD of PRESSURE (Pa^2/Hz), frequencies ff and df.

e   = e(:);
N   = numel(e);
df  = Fs / nfft;
Nf  = nfft/2 + 1;

% Build overlapping segments (50% overlap), detrend, Hanning window
numoverlap = nfft/2;
% number of complete segments with 50% overlap:
NS1 = floor(N/nfft);
if NS1 < 1, error('Not enough samples for nfft=%d', nfft); end
NS  = NS1*2 - 1;

zmat = zeros(nfft, NS);
zmat(:,1:NS1)         = reshape(e(1:nfft*NS1), [nfft, NS1]);
if NS > NS1
    zmat(:,NS1+1:NS)  = reshape(e(1+nfft/2 : nfft/2 + nfft*(NS1-1)), [nfft, NS1-1]);
end
zmat  = detrend(zmat,'linear');

% Hanning + bad-sample rejection (variance outliers), like your code
w    = hann(nfft,'periodic');
H    = w * ones(1,NS);
var2 = var(zmat); meanvar = mean(var2);
for it=1:4
    idx = find(var2 - 4*meanvar > 0);
    if isempty(idx), break; end
    zmat(:,idx) = 0;    % zero bad segments
    var2 = var(zmat); meanvar = mean(var2);
end
zmat  = H .* detrend(zmat,'constant');
wc2   = 1/mean(w.^2);                 % window correction
zspec = fft(zmat, nfft, 1)/nfft;
Szz   = mean(abs(zspec).^2, 2) * wc2 / df * 2;   % one-sided
PP    = Szz(2:Nf);                                % drop DC
ff    = (1:Nf-1).' * df;
end
