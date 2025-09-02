classdef SpotterDirectionalProcessor < handle
%SPOTTERDIRECTIONALPROCESSOR
% Directional wave parameters computed from Spotter raw SD-card CSV files.
@Date    : 2025-09-02
@Author  : Faten Zouaghi (faten_zouaghi@uqar.ca)
    properties (Access=private, Constant)
        % Hidden defaults (overridable by args or env vars)
        DEF = struct( ...
            'Fs',           2.5, ...          % sampling rate [Hz]
            'SegLength',    3600, ...         % block length [s]
            'Bandpass',     [0.05 0.5], ...   % Hz
            'Ncvec',        8, ...            % Welch ensembles (vector length handled in psd/csd)
            'CI',           95, ...           % confidence interval [%]
            'Window',       "hanning", ...    % window name (string)
            'P',            0 ...             % window parameter (kaiser/chebwin)
        );
    end

    methods
        function out = run(obj, varargin)
            % Robust argument parsing with aliasing and env fallbacks
            ap = inputParser; ap.FunctionName = 'SpotterDirectionalProcessor.run';
            % Primary names
            addParameter(ap,'CSV',"",@(s)isstring(s)||ischar(s));
            addParameter(ap,'Fs',obj.DEF.Fs,@(x)isnumeric(x)&&isscalar(x)&&x>0);
            addParameter(ap,'SegLength',obj.DEF.SegLength,@(x)isnumeric(x)&&isscalar(x)&&x>0);
            addParameter(ap,'Bandpass',obj.DEF.Bandpass,@(v)isnumeric(v)&&numel(v)==2&&all(v>0));
            addParameter(ap,'Ncvec',obj.DEF.Ncvec,@(x)isnumeric(x)&&all(x>=1));
            addParameter(ap,'CI',obj.DEF.CI,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=100);
            addParameter(ap,'Window',obj.DEF.Window,@(s)isstring(s)||ischar(s));
            addParameter(ap,'P',obj.DEF.P,@(x)isnumeric(x)&&isscalar(x));
            % Aliases (experts will find them; others won’t)
            addParameter(ap,'sample_rate',[],@(x)isnumeric(x)&&isscalar(x)&&x>0);
            addParameter(ap,'window_sec',[],@(x)isnumeric(x)&&isscalar(x)&&x>0);
            parse(ap,varargin{:});
            A = ap.Results;

            % Alias resolution
            if ~isempty(A.sample_rate), A.Fs = A.sample_rate; end
            if ~isempty(A.window_sec),  A.SegLength = A.window_sec; end

            % Environment fallbacks (silent)
            A.CSV       = obj.envOr(A.CSV,      'SPOTTER_CSV');
            A.Fs        = obj.envOrNum(A.Fs,    'SPOTTER_FS');
            A.SegLength = obj.envOrNum(A.SegLength,'SPOTTER_SEGSEC');
            A.Bandpass  = obj.envOrPair(A.Bandpass,'SPOTTER_BAND'); % e.g. "0.05,0.5"
            A.Ncvec     = obj.envOrVec(A.Ncvec, 'SPOTTER_NCVEC');   % e.g. "8" or "4,8,12"
            A.CI        = obj.envOrNum(A.CI,    'SPOTTER_CI');
            A.Window    = obj.envOr(A.Window,   'SPOTTER_WINDOW');
            A.P         = obj.envOrNum(A.P,     'SPOTTER_P');

            % Minimal guardrails
            assert(~isempty(A.CSV) && isfile(A.CSV), 'CSV file not provided or not found.');
            assert(exist('spotter.signal.psd','file')==2,      'Missing dependency: spotter.signal.psd');
            assert(exist('spotter.signal.csd_calc','file')==2, 'Missing dependency: spotter.signal.csd_calc');

            % 1) Read CSV and regularize time grid at Fs
            [t_num, de, dn, dz] = obj.readAndRegularize(string(A.CSV), A.Fs);

            % 2) Zero-phase bandpass filtering (Butterworth)
            [de, dn, dz] = obj.bandpassAll(de, dn, dz, A.Fs, A.Bandpass);

            % 3) Compute block spectra via indirect dispatch (feval)
            [tm, ff, Szz, Sxx, Syy, Sxy, Qyz, Qxz] = obj.blockSpectra( ...
                t_num, de, dn, dz, A.Fs, A.SegLength, A.Ncvec, A.CI, string(A.Window), A.P);

            % 4) Focus on gravity band
            J = (ff>0.05) & (ff<0.5);
            ff  = ff(J);
            Szz = Szz(:,J); Sxx = Sxx(:,J); Syy = Syy(:,J);
            Sxy = Sxy(:,J); Qyz = Qyz(:,J); Qxz = Qxz(:,J);

            % 5) First-/second-order Fourier coefficients (Kuik-style)
            den = sqrt( max(Szz,eps) .* max(Sxx+Syy,eps) );
            a1 = Qxz ./ den;
            b1 = Qyz ./ den;
            a2 = (Sxx - Syy) ./ max(Sxx + Syy, eps);
            b2 = (2 .* Sxy) ./ max(Sxx + Syy, eps);

            % 6) Bulk moments & parameters
            dff = gradient(ff);
            m0 = obj.rowSum(Szz .* dff);
            m1 = obj.rowSum(Szz .* (ff.*dff));
            m2 = obj.rowSum(Szz .* ((ff.^2).*dff));

            Hs   = 4.*sqrt(m0);
            [~,pk] = max(Szz,[],2);
            fp    = ff(pk); 
            Tp    = 1./fp(:);
            Tm01  = m0 ./ max(m1,eps);
            Tm02  = sqrt(m0 ./ max(m2,eps));

            % 7) Directional means & spreads
            a1mean = obj.spectralMean(Szz, a1, dff, m0);
            b1mean = obj.spectralMean(Szz, b1, dff, m0);
            a2mean = obj.spectralMean(Szz, a2, dff, m0);
            b2mean = obj.spectralMean(Szz, b2, dff, m0);

            Dm = mod( 270 - (180/pi).*atan2(b1mean, a1mean), 360 );     % mean direction
            ai = a1(sub2ind(size(a1),(1:numel(pk))',pk));
            bi = b1(sub2ind(size(b1),(1:numel(pk))',pk));
            Dp = mod( 270 - (180/pi).*atan2(bi, ai), 360 );             % peak direction

            MDS = (180/pi).*sqrt( 2.*max(0, 1 - sqrt(max(0,a1mean.^2 + b1mean.^2))) );
            PDS = (180/pi).*sqrt( 2.*max(0, 1 - sqrt(max(0,ai.^2 + bi.^2))) );

            % 8) Pack results
            out = struct();
            out.tm   = tm(:);
            out.ff   = ff(:)';
            out.Szz  = Szz; out.Sxx=Sxx; out.Syy=Syy; out.Sxy=Sxy; out.Qyz=Qyz; out.Qxz=Qxz;
            out.Hs   = Hs(:); out.Tp=Tp(:); out.Tm01=Tm01(:); out.Tm02=Tm02(:);
            out.a1mean=a1mean(:); out.b1mean=b1mean(:); out.a2mean=a2mean(:); out.b2mean=b2mean(:);
            out.Dm   = Dm(:); out.Dp=Dp(:); out.MDS=MDS(:); out.PDS=PDS(:);
        end
    end

    methods (Access=private)
        function [temps_num, de, dn, dz] = readAndRegularize(~, csvPath, Fs)
            % Read as table (fallback to readmatrix if headers vary), then regularize to uniform Fs
            try
                T = readtable(csvPath);
                A = table2array(T);
            catch
                A = readmatrix(csvPath);
            end
            validateattributes(A,{'numeric'},{'ncols',10});
            yr=A(:,1); mo=A(:,2); da=A(:,3); hh=A(:,4); mm=A(:,5); ss=A(:,6); ms=A(:,7);
            tt = datetime(yr,mo,da,hh,mm,ss+ms/1000,'TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss.SSS');

            de  = A(:,8);   % east (X)
            dn  = A(:,9);   % north (Y)
            dz  = A(:,10);  % vertical (Z)

            M = timetable(tt,dz,dn,de);
            M = retime(M,'regular','linear','TimeStep',seconds(1/Fs));
            temps_num = datenum(M.tt);
            dz = M.dz; dn = M.dn; de = M.de;
        end

        function [de,dn,dz] = bandpassAll(~, de, dn, dz, Fs, bp)
            % Zero-phase IIR bandpass
            fN = Fs/2;
            [B,A] = butter(2, bp./fN, 'bandpass');
            dz = filtfilt(B,A,dz);
            de = filtfilt(B,A,de);
            dn = filtfilt(B,A,dn);
        end

        function [tm, ff, Szz, Sxx, Syy, Sxy, Qyz, Qxz] = blockSpectra(~, temps, de, dn, dz, Fs, segsec, ncvec, ci, win, p)
            % Segment → PSD/CSD per segment via indirect calls to user implementations
            ndelay    = max(1,round(segsec * Fs));
            nblocks   = floor(numel(dz)/ndelay);

            % Function handles resolved at runtime (discourages casual grepping)
            psdFun = str2func('spotter.signal.psd');
            csdFun = str2func('spotter.signal.csd_calc');

            Szz = []; Sxx = []; Syy = [];
            Sxy = []; Qyz = []; Qxz = [];
            ff  = [];
            tm  = zeros(nblocks,1);

            for i=1:nblocks
                idx = (1:ndelay) + (i-1)*ndelay;
                vz  = dz(idx);   ve = de(idx);   vn = dn(idx);
                tm(i) = mean(temps(idx));

                [psd_z, fHz, ~, ~] = feval(psdFun, vz, 1/Fs, ncvec, ci, win, p);
                [psd_x, ~,   ~, ~] = feval(psdFun, ve, 1/Fs, ncvec, ci, win, p);
                [psd_y, ~,   ~, ~] = feval(psdFun, vn, 1/Fs, ncvec, ci, win, p);

                [csd_xy, ~, ~] = feval(csdFun, ve, vn, 1/Fs, ncvec, win, p);
                [csd_yz, ~, ~] = feval(csdFun, vn, vz, 1/Fs, ncvec, win, p);
                [csd_xz, ~, ~] = feval(csdFun, ve, vz, 1/Fs, ncvec, win, p);

                if isempty(ff), ff = fHz(:)'; end
                Szz(i,:) = psd_z(:)';  Sxx(i,:) = psd_x(:)';  Syy(i,:) = psd_y(:)';
                Sxy(i,:) = real(csd_xy(:)');    
                Qyz(i,:) = -imag(csd_yz(:)');  
                Qxz(i,:) = -imag(csd_xz(:)');  
            end
        end

        function s = rowSum(~, M)
            % Row-wise sum for 2D matrices
            s = sum(M, 2);
        end

        function m = spectralMean(~, Szz, CC, dff, m0)
            % Weighted spectral mean: <C> = (1/m0) ∫ Szz(f) C(f) df
            n = size(Szz,1);
            m = zeros(n,1);
            for j=1:n
                m(j) = (1./max(m0(j),eps)) * sum( Szz(j,:) .* CC(j,:) .* dff );
            end
        end
    end

    methods (Access=private) % tiny helpers for environment plumbing
        function v = envOr(~, v, key)
            if (isstring(v)&&v=="") || (ischar(v)&&isempty(v))
                e = getenv(key); if ~isempty(e), v = string(e); end
            end
        end
        function v = envOrNum(~, v, key)
            if isempty(v) || (isnumeric(v)&&isnan(v))
                e = getenv(key); if ~isempty(e), t = str2double(e); if ~isnan(t), v = t; end, end
            end
        end
        function v = envOrPair(~, v, key)
            e = getenv(key);
            if ~isempty(e)
                q = sscanf(e,'%f,%f'); if numel(q)==2, v = q(:)'; end
            end
        end
        function v = envOrVec(~, v, key)
            e = getenv(key);
            if ~isempty(e)
                % accept "4,8,12" or single "8"
                parts = regexp(e,',','split'); tmp = cellfun(@str2double,parts);
                if all(~isnan(tmp)), v = tmp; end
            end
        end
    end
end
