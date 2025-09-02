@Date    : 2025-09-02
@Author  : Faten Zouaghi (faten_zouaghi@uqar.ca)

classdef SpotterDirectionalProcessor
% SpotterDirectionalProcessor
% Compute directional wave parameters from Spotter raw SD-card CSV files.
%

    properties (Access=private, Constant)
        % Required keys: no embedded numeric defaults. User must supply them.
        REQ_KEYS = {'Fs','SegLength','Bandpass','Ncvec','CI','Window','P'};
    end

    methods
        function out = run(obj, varargin)
            % ---- Parse without numeric defaults (forces expert input) ----
            ap = inputParser; ap.FunctionName = 'SpotterDirectionalProcessor.run';
            addParameter(ap,'CSV',"",@(s)isstring(s)||ischar(s));
            addParameter(ap,'Fs',[],@(x)isnumeric(x)&&isscalar(x)&&x>0);
            addParameter(ap,'SegLength',[],@(x)isnumeric(x)&&isscalar(x)&&x>0);
            addParameter(ap,'Bandpass',[],@(v)isnumeric(v)&&numel(v)==2&&all(v>0));
            addParameter(ap,'Ncvec',[],@(x)isnumeric(x)&&all(x>=1));
            addParameter(ap,'CI',[],@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=100);
            addParameter(ap,'Window',"",@(s)isstring(s)||ischar(s));
            addParameter(ap,'P',[],@(x)isnumeric(x)&&isscalar(x));
            % hidden aliases (also empty by default)
            addParameter(ap,'sample_rate',[],@(x)isnumeric(x)&&isscalar(x)&&x>0);
            addParameter(ap,'window_sec',[],@(x)isnumeric(x)&&isscalar(x)&&x>0);
            parse(ap,varargin{:});
            A = ap.Results;

            % Alias resolution
            if ~isempty(A.sample_rate), A.Fs = A.sample_rate; end
            if ~isempty(A.window_sec),  A.SegLength = A.window_sec; end

            % ---- Expert-only fallbacks: ENV → MATLAB preferences ----
            A.CSV       = obj.i_envOr(A.CSV,      'SPOTTER_CSV');
            A.Fs        = obj.i_envOrNum(A.Fs,    'SPOTTER_FS');
            A.SegLength = obj.i_envOrNum(A.SegLength,'SPOTTER_SEGSEC');
            A.Bandpass  = obj.i_envOrPair(A.Bandpass,'SPOTTER_BAND');   % e.g. "0.05,0.5"
            A.Ncvec     = obj.i_envOrVec(A.Ncvec, 'SPOTTER_NCVEC');     % e.g. "8" or "4,8,12"
            A.CI        = obj.i_envOrNum(A.CI,    'SPOTTER_CI');
            A.Window    = obj.i_envOr(A.Window,   'SPOTTER_WINDOW');
            A.P         = obj.i_envOrNum(A.P,     'SPOTTER_P');

            A.Fs        = obj.i_prefOr('spotter','Fs',A.Fs);
            A.SegLength = obj.i_prefOr('spotter','SegLength',A.SegLength);
            A.Bandpass  = obj.i_prefOr('spotter','Bandpass',A.Bandpass);
            A.Ncvec     = obj.i_prefOr('spotter','Ncvec',A.Ncvec);
            A.CI        = obj.i_prefOr('spotter','CI',A.CI);
            A.Window    = obj.i_prefOr('spotter','Window',A.Window);
            A.P         = obj.i_prefOr('spotter','P',A.P);

            % ---- Guards (blocks non-expert usage) ----
            assert(~isempty(A.CSV) && isfile(A.CSV), ...
                'CSV file not provided or not found (set ''CSV'' or SPOTTER_CSV).');
            missing = {};
            for k = 1:numel(obj.REQ_KEYS)
                key = obj.REQ_KEYS{k}; val = A.(key);
                if (isnumeric(val) && (isempty(val) || any(isnan(val)))) || ...
                   (isstring(val) && strlength(val)==0) || (ischar(val) && isempty(val))
                    missing{end+1} = key; %#ok<AGROW>
                end
            end
            assert(isempty(missing), 'Missing required settings: %s', strjoin(missing, ', '));

            % Required dependencies (kept under package for obscurity)
            assert(exist('spotter.signal.psd','file')==2, 'Missing dependency: src/+spotter/+signal/psd.m');
            assert(exist('spotter.signal.csd','file')==2, 'Missing dependency: src/+spotter/+signal/csd.m');

            % ---- 1) Read CSV & regularize to Fs ----
            [t_num, de, dn, dz] = obj.readAndRegularize(string(A.CSV), A.Fs);

            % ---- 2) Band-pass filter on all components ----
            [de, dn, dz] = obj.bandpassAll(de, dn, dz, A.Fs, A.Bandpass);

            % ---- 3) Block-wise spectra & cross-spectra ----
            [tm, ff, Szz, Sxx, Syy, Sxy, Qyz, Qxz] = obj.blockSpectra( ...
                t_num, de, dn, dz, A.Fs, A.SegLength, A.Ncvec, A.CI, string(A.Window), A.P);

            % ---- 4)  band selection (0.05–0.5 Hz) ----
            J   = (ff > 0.05) & (ff < 0.5);
            ff  = ff(J);
            Szz = Szz(:,J); Sxx = Sxx(:,J); Syy = Syy(:,J);
            Sxy = Sxy(:,J); Qyz = Qyz(:,J); Qxz = Qxz(:,J);

            % ---- 5) Fourier coefficients ----
            den = sqrt( Szz .* (Sxx + Syy) );
            a1 = Qxz ./ den;                 
            b1 = Qyz ./ den;
            a2 = (Sxx - Syy) ./ (Sxx + Syy);
            b2 = (2 .* Sxy)  ./ (Sxx + Syy);

            % ---- 6) Bulk parameters
            dff = gradient(ff);
            m0  = obj.rowSum(Szz .* dff);
            m1  = obj.rowSum(Szz .* (ff.*dff));
            m2  = obj.rowSum(Szz .* ((ff.^2).*dff));

            Hs   = 4.*sqrt(m0);
            [~,pk] = max(Szz,[],2);
            fp    = ff(pk);  Tp   = 1./fp(:);
            Tm01  = m0 ./ max(m1,eps);
            Tm02  = sqrt(m0 ./ max(m2,eps));

        
            a1mean = obj.spectralMean(Szz, a1, dff, m0);
            b1mean = obj.spectralMean(Szz, b1, dff, m0);
            a2mean = obj.spectralMean(Szz, a2, dff, m0);
            b2mean = obj.spectralMean(Szz, b2, dff, m0);

     
            Dm  = mod( 270 - (180/pi).*atan2(b1mean, a1mean), 360 );
            ai  = a1(sub2ind(size(a1),(1:numel(pk))',pk));
            bi  = b1(sub2ind(size(b1),(1:numel(pk))',pk));
            Dp  = mod( 270 - (180/pi).*atan2(bi, ai), 360 );

            MDS = mod( (180/pi).*sqrt( 2.*(1 - sqrt(a1mean.^2 + b1mean.^2)) ), 360 );
            PDS = mod( (180/pi).*sqrt( 2.*(1 - sqrt(ai.^2 + bi.^2)) ), 360 );

            % ---- Pack output ----
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
        function [t_num, de, dn, dz] = readAndRegularize(~, csvPath, Fs)
            % Read raw Spotter displacement CSV and regularize at Fs Hz.
            T = readtable(csvPath);
            A = table2array(T);
            yr=A(:,1); mo=A(:,2); da=A(:,3); hh=A(:,4); mm=A(:,5); ss=A(:,6); ms=A(:,7);
            tt = datetime(yr,mo,da,hh,mm,ss+ms/1000,'TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss.SSS');

            de = A(:,8);   % East (X)
            dn = A(:,9);   % North (Y)
            dz = A(:,10);  % Vertical (Z)

            M  = timetable(tt,dz,dn,de);
            M  = retime(M,'regular','linear','TimeStep',seconds(1/Fs));
            t_num = datenum(M.tt);
            dz = M.dz; dn = M.dn; de = M.de;
        end

        function [de,dn,dz] = bandpassAll(~, de, dn, dz, Fs, bp)
            % Apply 2nd-order Butterworth band-pass 
            fN = Fs/2;
            [B,A] = butter(2, bp./fN, 'bandpass');
            dz = filtfilt(B,A,dz);
            de = filtfilt(B,A,de);
            dn = filtfilt(B,A,dn);
        end

        function [tm, ff, Szz, Sxx, Syy, Sxy, Qyz, Qxz] = blockSpectra(~, t, de, dn, dz, Fs, segsec, ncvec, ci, win, p)
            % Split into blocks and compute PSD/CSD per block using package functions.
            ndelay    = round(segsec * Fs);
            nblocks   = floor(numel(dz)/ndelay);

            Szz = []; Sxx = []; Syy = [];
            Sxy = []; Qyz = []; Qxz = [];
            ff  = []; tm  = zeros(nblocks,1);

            for i=1:nblocks
                idx  = (1:ndelay) + (i-1)*ndelay;
                vz   = dz(idx);
                ve   = de(idx);
                vn   = dn(idx);
                tm(i)= mean(t(idx));

                dt = 1/Fs;
                [psd_z, fHz, ~, ~] = spotter.signal.psd(vz, dt, ncvec, ci, win, p);
                [psd_x, ~,   ~, ~] = spotter.signal.psd(ve, dt, ncvec, ci, win, p);
                [psd_y, ~,   ~, ~] = spotter.signal.psd(vn, dt, ncvec, ci, win, p);

                [csd_xy, ~, ~] = spotter.signal.csd(ve, vn, dt, ncvec, win, p);
                [csd_yz, ~, ~] = spotter.signal.csd(vn, vz, dt, ncvec, win, p);
                [csd_xz, ~, ~] = spotter.signal.csd(ve, vz, dt, ncvec, win, p);

                if isempty(ff), ff = fHz(:)'; end
                Szz(i,:) = psd_z(:)';  Sxx(i,:) = psd_x(:)';  Syy(i,:) = psd_y(:)';
                Sxy(i,:) = real(csd_xy(:)');    
                Qyz(i,:) = -imag(csd_yz(:)');  
                Qxz(i,:) = -imag(csd_xz(:)');  
            end
        end

        function s = rowSum(~, M)
            % Row-wise sum for 2D matrices.
            s = sum(M, 2);
        end

        function m = spectralMean(~, Szz, C, dff, m0)
            % Weighted spectral mean: sum(Szz .* C .* dff) / m0 (per row).
            n = size(Szz,1);
            m = zeros(n,1);
            for j=1:n
                m(j) = (1./m0(j)) * sum( Szz(j,:) .* C(j,:) .* dff );
            end
        end


        function v = i_envOr(~, v, key)
            if (isstring(v)&&v=="") || (ischar(v)&&isempty(v))
                e = getenv(key); if ~isempty(e), v = string(e); end
            end
        end

        function v = i_envOrNum(~, v, key)
            if isempty(v) || (isnumeric(v)&&any(isnan(v)))
                e = getenv(key);
                if ~isempty(e)
                    t = str2double(e);
                    if ~isnan(t), v = t; end
                end
            end
        end

        function v = i_envOrPair(~, v, key)
            e = getenv(key);
            if ~isempty(e)
                q = sscanf(e,'%f,%f');
                if numel(q)==2, v = q(:)'; end
            end
        end

        function v = i_envOrVec(~, v, key)
            e = getenv(key);
            if ~isempty(e)
                parts = regexp(e,',','split'); tmp = cellfun(@str2double,parts);
                if all(~isnan(tmp)), v = tmp; end
            end
        end

        function v = i_prefOr(~, ns, name, v)
            if (isnumeric(v) && (isempty(v) || any(isnan(v)))) || ...
               (ischar(v) && isempty(v)) || (isstring(v) && strlength(v)==0)
                try
                    if ispref(ns, name), v = getpref(ns, name); end
                catch, end
            end
        end
    end
end
