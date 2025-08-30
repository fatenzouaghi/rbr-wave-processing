% Code to obtain breaking statistics from Richard Manasseh's data
function [ff, df, PP, fp, Hs] = spectra (e, Fs, nfft, d)

z   = e;
nmes = length(e);

% We remove some trends
maxe = max(z);

NSbad = 0;
%Fs=25.;                        % sampling frequency in Herz
%nmes=ind2(irec);
% nfft=256.;                    % nfft specifies the FFT length that csd uses.

df       = Fs/nfft;               % Frequential resolution 
numoverlap = nfft/2;              % number of samples by which the sections overlap
Nf       = nfft/2 + 1;
Nfo      = (Nf - 1)/2;

See = zeros(Nf);                  % Power spectrum 
Ef  = zeros(Nf - 1);              % Power spectrum density

hanning = transpose(0.5*(1 - cos(2*pi*linspace(0,nfft-1,nfft)/(nfft-1))));
NS1     = floor(nmes/nfft);
NS      = NS1*2 - 1;
zmat    = zeros(nfft,NS);

zmat(1:nfft,1:NS1)      = reshape(z(1:nfft*NS1), [nfft,NS1]);
zmat(1:nfft,NS1+1:NS)   = reshape(z(1+nfft/2 : nfft/2 + nfft*(NS1-1)), [nfft,NS1-1]);

zmat2   = detrend(zmat,'linear');  % remove linear trend from each column

% PLOT OF DETRENDED TIME SERIES
Hs1 = 4.*sqrt(mean(var(zmat2)));
Eh  = ones(1,NS);
H   = hanning*Eh;

NSbad = 0;
nstep = 4;  % checks for bad samples
while (nstep > 0)
    nstep = nstep - 1;
    var2 = var(zmat2);
    meanvar = mean(var2);
    index = find(var2 - 4*meanvar > 0);
    taille = size(index);
    NSbad = NSbad + taille(2);
    if (NSbad == 0)
        nstep = 0;
    end
    zmat2(:,index) = 0;
end

zmat2 = H.*detrend(zmat2,'constant');   % remove mean and apply window
Hs2   = 4.*sqrt(mean(var(zmat2)));
wc2   = 1/mean(hanning.^2);             % window correction factor
fac   = NS/(NS-NSbad);
zspec = fft(zmat2,nfft,1)/nfft;
Szzmat = (abs(zspec).^2)*wc2/df;
Czz    = mean((abs(zspec).^2),2)*wc2/df*2;
freq1  = linspace(0, df*(nfft-1), nfft);
%plot(freq1,transpose(Seemat(1,:)))

% converts two sided spectral density to one sided spectral density
Szzmat = Szzmat.*2; 
depth  = d;
Szz    = mean(Szzmat,2)*fac;
freq   = linspace(df, df*(Nf-1), Nf-1); 
Ef     = transpose(Szz(2:Nf));
%Ef(1) = Ef(2); % removes low frequency variance

Hs = 4.*sqrt(sum(Ef(:))*df);
[Emax, ifp] = max(Ef(1:(Nf-1)/2));

PP  = Ef(2:end);
ff  = freq(2:end);
dff = diff(ff);
df  = dff(1);

[YY, I] = max(PP(3:end));
fp = ff(I+2); % peak frequency

end

