function [psd,freq,errlo,errhi] = psd_calc(x,dt,ncvec,ci,win_name,p)
% [psd,freq,errlo,errhi] = psd_calc(x,dt,ncvec,ci,win_name,p)
%
% Estimates the power spectral density of a regularly spaced time series,
% using Welch's averaged periodogram method (ensemble-averaging)
% with variable number of ensembles to be specified.
%
% INPUT: x = data record (vector) [xunit]
%	 dt = time step (scalar) [timeunit]
%	 ncvec = list of number of ensembles to use for averaging (vector).
%                The logarithmic frequency axis is divided into length(ncvec) sections.
%                To compute the psd for the frequencies in section i, the time series is 
%                divided into ncvec(i) half-overlapping segments.
%        ci = confidence interval [%] for the error estimates.
%        win_name = name of window to be applied to each segment of the time series
%                   before the psd is computed ('bartlett','blackman','boxcar','chebwin',
%                   'hamming','hann','hanning','kaiser','triang').
%        p = window parameter (only useful for windows kaiser and chebwin, but put
%            something even if you use another window otherwise the program will not run).
% 
% OUTPUT: psd = power spectral density (column vector) [xunit^2 / cycles-per-timeunit]
%         freq = frequency (column vector) [cycles-per-timeunit]
%	  errlo = 95% confidence lower limit (column vector)
%	  errhi = 95% confidence upper limit (column vector)
%
% NOTE: you multiply the signal or some constant by [errlo,errhi] to put
%	them on the loglog plot

% psd_calc.m 03/25/2005 Cedric Chavanne, adapted from
% psd_calc.m 1/6/2000 and psd_plot.m 3/19/2002 Parker MacCready
% 12/19/2013 modified the detrend step so that the mean is also removed

nn = length(ncvec);
% set frequency limits
fhi = 1/(2*dt); % the highest (Nyquist) frequency
% It is a little trickier to find the low frequency limit
% because it depends on the number of blocks used in the lowest
% frequency range. We do it exactly as in "psd_calc.m"
N = length(x);  % initial record length
nc = ncvec(1); % the number of blocks in the lowest frequency range
nr2 = floor(N/(nc+1));  % the maximum length of each half block
nr = 2*nr2;
flo = 1/(nr*dt)-eps;    % the lowest frequency of the averaged signal
% we substract "eps" because otherwise the "mask" used below seems to
% miss the lowest frequency point.
%
% and make the frequency break points
break_fr = logspace(log10(flo),log10(fhi),nn+1);

% set up the vectors to hold the results of all the
% different block averages
psd_all = [];
freq_all = [];
elo_all = [];
ehi_all = [];
dof_all = [];

% now do the psd calculations and put the results into the
% parts of the result vectors.  We use "break_fr" to divvy up
% the results evenly in logscale frequency space
for ii = 1:nn
  nc = ncvec(ii);
  % reshape the signal to have "nc" columns, each having 50% overlap
  nr2 = floor(N/(nc+1));	% the maximum length of each half block
  nr = 2*nr2;
  % make the reshaped signal matrix "xx"
  xx = [];
  x = x(:);
  for mm = 1:nc
	nstart = 1 + (mm-1)*nr2;
	nend = nstart + nr - 1;
	xx = [xx, x(nstart:nend)];
  end
  % and find the new total record length
  Nii = nr2 * (nc+1);
  %
  % compute the psd of each column using a Hanning window
  %	(should work for arbitrary window shape)
%  win = hanning(nr);
  if strcmp(win_name,'kaiser')==1 | strcmp(win_name,'chebwin')==1
    eval(['win = ',win_name,'(nr,p);']);
  else
    eval(['win = ',win_name,'(nr);']);
  end
  win_mat = win * ones(1,nc);
%  XX = fft(win_mat.*detrend(xx));
  XX = fft(win_mat.*detrend(detrend(xx,'constant')));
%  psd = 2 * dt * XX(2:nr2+1,:) .* conj(XX(2:nr2+1,:)) / sum(win.^2);
  psd = dt * abs(XX).^2 / sum(win.^2);
  % and now average aross columns (i.e. average each frequency band)
  psd = mean(psd,2);
%  psd = real(psd);	% the imaginary part is all zero
  % and make the frequency vector
%  freq = [1:nr2]' / (nr*dt);
  freq = [0:nr2 , -nr2+1:-1]' / (nr*dt);
  % and calculate the 95% confidence limits
  dof = (Nii/nr2) * (nr+1) / sum(win.^2);
%  [elo,ehi] = chisqp(dof,95); % "chisqp.m" is a Dewey routine
  alpha = 1-ci/100;
  elo = dof/chi2inv(1-alpha/2,dof);
  ehi = dof/chi2inv(alpha/2,dof);
   elo = elo*ones(size(freq));
   ehi = ehi*ones(size(freq));
   dof = dof*ones(size(freq));
   this_flo = break_fr(ii)-eps;
   this_fhi = break_fr(ii+1);
%   mask = (freq>=this_flo & freq<this_fhi);
   maskp = (freq>=this_flo & freq<this_fhi);
   maskn = (freq<=-this_flo & freq>-this_fhi);
%   psd_all = [psd_all; psd(mask)];
%   freq_all = [freq_all; freq(mask)];
%   elo_all = [elo_all; elo(mask)];
%   ehi_all = [ehi_all; ehi(mask)];
   psd_all = [psd(maskn);psd_all; psd(maskp)];
   freq_all = [freq(maskn);freq_all; freq(maskp)];
   elo_all = [elo(maskn);elo_all; elo(maskp)];
   ehi_all = [ehi(maskn);ehi_all; ehi(maskp)];
   dof_all = [dof(maskn);dof_all; dof(maskp)];

%output
psd=psd_all;
freq=freq_all;
errlo=elo_all;
errhi=ehi_all;
dof=dof_all;
end
