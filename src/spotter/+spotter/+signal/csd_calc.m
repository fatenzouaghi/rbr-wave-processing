function [csd,freq,dof] = csd_calc(x,y,dt,ncvec,win_name,p)

% [csd,freq,dof] = csd_calc(x,y,dt,ncvec,win_name,p)
%
% Estimates the cross spectral density of two regularly spaced time series,
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
csd_all = [];
freq_all = [];
dof_all = [];

% now do the psd calculations and put the results into the
% parts of the result vectors.  We use "break_fr" to divvy up
% the results evenly in logscale frequency space
for ii = 1:nn
  nc = ncvec(ii);
  % reshape the signal to have "nc" columns, each having 50% overlap
  nr2 = floor(N/(nc+1));	% the maximum length of each half block
  nr = 2*nr2;
  % make the reshaped signal matrices "xx" and "yy"
  xx = []; yy = [];
  x = x(:); y = y(:);
  for mm = 1:nc
	nstart = 1 + (mm-1)*nr2;
	nend = nstart + nr - 1;
	xx = [xx, x(nstart:nend)];
	yy = [yy, y(nstart:nend)];
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
  XX = fft(win_mat.*detrend(detrend(xx,'constant')));
  YY = fft(win_mat.*detrend(detrend(yy,'constant')));
  csd = dt * conj(XX) .* YY / sum(win.^2);
  % and now average aross columns (i.e. average each frequency band)
  csd = mean(csd,2);
  % and make the frequency vector
%  freq = [1:nr2]' / (nr*dt);
  freq = [0:nr2 , -nr2+1:-1]' / (nr*dt);
  % and calculate the degrees of freedom:
  dof = (Nii/nr2) * (nr+1) / sum(win.^2);
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
   csd_all = [csd(maskn);csd_all; csd(maskp)];
   freq_all = [freq(maskn);freq_all; freq(maskp)];
   dof_all = [dof(maskn);dof_all; dof(maskp)];
end
%output
csd=csd_all;
freq=freq_all;
dof=dof_all;

